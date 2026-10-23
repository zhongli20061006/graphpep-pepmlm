#!/usr/bin/env python3
"""Score each real complex with its own receptor. Raw GraphPep score: lower better.

Different sequences cannot share a multi-MODEL topology or ESM feature. Optional
manifest columns: complex_path, model_id, receptor_sequence. Confidence columns
are preserved. Empty interaction graphs stay unscored.
"""
from __future__ import annotations
import argparse
import csv
import math
from pathlib import Path
import re
import subprocess
import sys
import time
import numpy as np
import pandas as pd

ROOT = Path("/home/g/pepmlm_pipeline")
GP = Path("/home/g/GraphPep_v1.3")
BLUE_THRESHOLD = 0.70
THREE_TO_ONE = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))


def read_models(path: str | Path) -> list[list[str]]:
    """Read PDB models without mixing atom records."""
    models, current = [], []
    for line in Path(path).read_text().splitlines(keepends=True):
        if line.startswith("MODEL"):
            if current:
                models.append(current)
            current = []
        elif line.startswith("ENDMDL"):
            models.append(current)
            current = []
        elif line.startswith(("ATOM", "HETATM", "TER")):
            current.append(line)
    if current:
        models.append(current)
    return models


def split_chain(lines: list[str], chain: str) -> tuple[list[str], str, dict[int, str]]:
    """Normalize residue IDs while preserving heavy atoms and original IDs."""
    selected, seen_atoms, residue_ids, sequence, mapping = [], set(), {}, [], {}
    for line in lines:
        if not line.startswith("ATOM") or len(line) < 54 or line[21] != chain:
            continue
        name = line[12:16].strip()
        element = line[76:78].strip() if len(line) >= 78 else ""
        if element in {"H", "D"} or (not element and name.lstrip("0123456789").startswith("H")):
            continue
        if line[16] not in {" ", "A"}:
            continue
        residue = (line[22:26].strip(), line[26])
        key = (*residue, name)
        if key in seen_atoms:
            continue
        seen_atoms.add(key)
        resname = line[17:20].strip()
        if resname not in THREE_TO_ONE:
            raise ValueError(f"unsupported residue {resname} in chain {chain}")
        if residue not in residue_ids:
            index = len(residue_ids) + 1
            residue_ids[residue] = index
            sequence.append(THREE_TO_ONE[resname])
            mapping[index] = f"{residue[0]}{residue[1].strip()}:{chain}"
        index = residue_ids[residue]
        normalized = line[:16] + " " + line[17:22] + f"{index:4d} " + line[27:]
        selected.append(normalized.rstrip("\r\n") + "\n")
    if not selected:
        raise ValueError(f"no standard protein atoms in chain {chain}")
    return selected, "".join(sequence), mapping


def build_inputs(manifest: pd.DataFrame, pdb_dir: str) -> list[dict]:
    """Produce a separate coordinate-preserving pair for each manifest row."""
    base, pairs = Path(pdb_dir), []
    for i, row in enumerate(manifest.to_dict("records")):
        ident = str(row["peptide_id"])
        value = row.get("complex_path")
        cf = Path(str(value)) if pd.notna(value) and value else base / f"{ident}_complex.pdb"
        if not cf.is_file():
            raise FileNotFoundError(f"missing complex: {cf}")
        models = read_models(cf)
        if len(models) != 1:
            raise ValueError(f"{cf}: expected one model; import structures separately")
        rec, rec_seq, rec_map = split_chain(models[0], "A")
        pep, pep_seq, pep_map = split_chain(models[0], "B")
        if pep_seq != str(row["peptide_sequence"]):
            raise ValueError(f"{cf}: peptide sequence {pep_seq} != manifest sequence {row['peptide_sequence']}")
        expected_rec = row.get("receptor_sequence")
        if pd.notna(expected_rec) and expected_rec and rec_seq != expected_rec:
            raise ValueError(f"{cf}: receptor sequence does not match manifest")
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{i:04d}_{ident}_{row.get('model_id', 'model1')}")
        folder = base / "graphpep_inputs" / name
        folder.mkdir(parents=True, exist_ok=True)
        rec_path, pep_path = folder / "protein.pdb", folder / "decoys.pdb"
        rec_path.write_text("".join(rec) + "TER\nEND\n")
        pep_path.write_text("MODEL        1\n" + "".join(pep) + "TER\nENDMDL\nEND\n")
        with (folder / "residue_mapping.csv").open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["graphpep_residue", "original_residue"])
            for ch, mapping in [("A", rec_map), ("B", pep_map)]:
                writer.writerows((f"{n}:{ch}", original) for n, original in mapping.items())
        pairs.append({"row": row, "folder": folder, "receptor_path": rec_path,
                      "peptide_path": pep_path, "receptor_map": rec_map, "peptide_map": pep_map})
    return pairs


def contact_count(pair: dict, cutoff: float = 6.) -> int:
    """Count heavy-atom residue contacts using the model's distance definition."""
    def residues(path):
        groups = {}
        for line in Path(path).read_text().splitlines():
            if line.startswith("ATOM"):
                groups.setdefault(int(line[22:26]), []).append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])])
        return {index: np.asarray(points) for index, points in groups.items()}
    rec, pep = residues(pair["receptor_path"]), residues(pair["peptide_path"])
    return sum(bool(np.any(np.sum((a[:, None, :] - b[None, :, :]) ** 2, axis=2) <= cutoff ** 2))
               for a in rec.values() for b in pep.values())


def map_feedback(each: pd.DataFrame, pair: dict) -> pd.DataFrame:
    """Restore original IDs and validate contact probabilities."""
    required = {"DECOY_ID", "RESID_PEP", "RESID_PRO", "CONFIDENT_SCORES"}
    if not required.issubset(each.columns):
        raise ValueError(f"invalid feedback columns: {list(each.columns)}")
    each = each.copy()
    each["CONFIDENT_SCORES"] = pd.to_numeric(each["CONFIDENT_SCORES"], errors="raise")
    if not each["CONFIDENT_SCORES"].between(0, 1).all():
        raise ValueError("contact probability outside [0, 1]")
    for field, mapping, ch in [("RESID_PEP", pair["peptide_map"], "B"), ("RESID_PRO", pair["receptor_map"], "A")]:
        each[f"GRAPHPEP_{field}"] = each[field]
        original = []
        for value in each[field].astype(str):
            num, chain = value.split(":")
            if chain != ch:
                raise ValueError(f"unexpected chain in {field}: {value}")
            original.append(mapping[int(num)])
        each[field] = original
    each["peptide_position"] = each["GRAPHPEP_RESID_PEP"].astype(str).str.split(":").str[0].astype(int)
    each["peptide_id"] = pair["row"]["peptide_id"]
    each["model_id"] = pair["row"].get("model_id", "model1")
    return each


def score_pair(pair: dict, graphpep_dir: Path, patched_resmap: Path, threshold: float) -> dict:
    """Invoke the installed program in a separate folder, retaining raw results."""
    output = dict(pair["row"])
    output.update(binding_score=float("nan"), red_residue_ratio=float("nan"),
                  low_conf_contact_pair_ratio=float("nan"), n_contacts=0,
                  min_conf=float("nan"), mean_conf=float("nan"), graphpep_status="pending")
    count = contact_count(pair)
    output.update(n_geometry_contact_pairs=count, graphpep_workdir=str(pair["folder"]))
    if count == 0:
        output["graphpep_status"] = "no_contact"
        return output
    # Upstream shell errors can still end with exit 0. Never consume stale outputs.
    for name in ["score.csv", "score_each.csv", "contact_feedback.csv", "residue_feedback.csv"]:
        (pair["folder"] / name).unlink(missing_ok=True)
    if not patched_resmap.is_file():
        raise FileNotFoundError(f"missing pipeline residue mapper: {patched_resmap}")
    shell = (graphpep_dir / "GraphPep.sh").read_text().replace("$GraphPep_root/bin/resmap.py", str(patched_resmap))
    (pair["folder"] / "graphpep_patched.sh").write_text(shell)
    t0 = time.perf_counter()
    proc = subprocess.run(["bash", "graphpep_patched.sh", "protein.pdb", "decoys.pdb", "-out", "score.csv",
                           "-each", "-outeach", "score_each.csv"], cwd=pair["folder"], capture_output=True, text=True)
    output["graphpep_seconds"] = time.perf_counter() - t0
    (pair["folder"] / "run.log").write_text(proc.stdout + "\nSTDERR\n" + proc.stderr)
    if proc.returncode:
        raise RuntimeError(f"GraphPep exit {proc.returncode}; see {pair['folder'] / 'run.log'}")
    score = pd.read_csv(pair["folder"] / "score.csv")
    if len(score) != 1 or int(score.iloc[0]["DECOY_ID"]) != 1:
        raise ValueError("expected one explicitly identified score")
    value = float(score.iloc[0]["PRED_SCORE"])
    if not math.isfinite(value):
        raise ValueError("score is not finite")
    each = map_feedback(pd.read_csv(pair["folder"] / "score_each.csv"), pair)
    if len(each) != count:
        raise ValueError(f"contact count mismatch: geometry {count}, GraphPep {len(each)}")
    each_path = pair["folder"] / "contact_feedback.csv"
    each.to_csv(each_path, index=False)
    positions = each.groupby("peptide_position")["CONFIDENT_SCORES"].agg(["count", "mean", "min", "max"])
    positions = positions.reindex(range(1, len(str(pair["row"]["peptide_sequence"])) + 1))
    positions["count"] = positions["count"].fillna(0).astype(int)
    positions.to_csv(pair["folder"] / "residue_feedback.csv")
    conf = each["CONFIDENT_SCORES"]
    ratio = float((conf < threshold).mean())
    output.update(binding_score=value, red_residue_ratio=ratio, low_conf_contact_pair_ratio=ratio,
                  n_contacts=len(each), min_conf=float(conf.min()), mean_conf=float(conf.mean()),
                  contact_feedback_path=str(each_path), residue_feedback_path=str(pair["folder"] / "residue_feedback.csv"),
                  graphpep_status="ok")
    return output


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lists", default=str(ROOT / "output/eval_lists.csv"))
    ap.add_argument("--pdb-dir", default=str(ROOT / "output/pdbs"))
    ap.add_argument("--scores", default=str(ROOT / "output/graphpep_scores.csv"))
    ap.add_argument("--graphpep-dir", type=Path, default=GP)
    ap.add_argument("--patched-resmap", type=Path, default=ROOT / "pipeline/graphpep_patches/resmap_compat.py")
    ap.add_argument("--blue-threshold", type=float, default=BLUE_THRESHOLD)
    args = ap.parse_args()
    if not args.patched_resmap.is_file():
        subprocess.run([sys.executable, str(ROOT / "pipeline/make_resmap_compat.py")], check=True)
    manifest = pd.read_csv(args.lists)
    if manifest.empty:
        raise ValueError("no structures in manifest; cloud input is pending")
    pairs, results = build_inputs(manifest, args.pdb_dir), []
    out = Path(args.scores)
    out.parent.mkdir(parents=True, exist_ok=True)
    for pair in pairs:
        print(f"scoring {pair['row']['peptide_id']} / {pair['row'].get('model_id', 'model1')}", flush=True)
        result = score_pair(pair, args.graphpep_dir, args.patched_resmap, args.blue_threshold)
        results.append(result)
        pd.DataFrame(results).to_csv(out, index=False)
        print(f"  {result['graphpep_status']}: {result['binding_score']}", flush=True)
    print(f"wrote {len(results)} model rows: {out}; lower raw score is better")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
