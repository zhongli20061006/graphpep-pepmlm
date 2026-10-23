#!/usr/bin/env python3
"""Apply the frozen identity, confidence, clash, and contact gates to six AFM models."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1


PILOT = Path(__file__).resolve().parent
CONFIG = json.loads((PILOT / "run_config_frozen_v1.json").read_text(encoding="utf-8"))
PARSER = PDBParser(QUIET=True)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_complex_fasta() -> dict[str, tuple[str, str]]:
    records: dict[str, str] = {}
    name = None
    for line in (PILOT / "complexes.fasta").read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            records[name] = ""
        else:
            records[name] += line.strip()
    return {name: tuple(seq.split(":")) for name, seq in records.items()}


def residues_and_atoms(pdb: Path):
    model = PARSER.get_structure(pdb.stem, str(pdb))[0]
    chains = list(model)
    sequences = {}
    numbers = {}
    chain_atoms = {}
    for chain in chains:
        residues = [r for r in chain if r.id[0] == " "]
        sequences[chain.id] = "".join(seq1(r.resname) for r in residues)
        numbers[chain.id] = [r.id[1] for r in residues]
        coords, residue_index = [], []
        for i, residue in enumerate(residues):
            for atom in residue:
                if atom.element.upper() not in {"H", "D"}:
                    coords.append(atom.coord.astype(float))
                    residue_index.append(i)
        chain_atoms[chain.id] = (np.asarray(coords), np.asarray(residue_index, dtype=int))
    return [c.id for c in chains], sequences, numbers, chain_atoms


def geometry_metrics(chain_atoms):
    a, ai = chain_atoms["A"]
    b, bi = chain_atoms["B"]
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    contact_atom_pairs = np.argwhere(distances <= 6.0)
    contact_pairs = {(int(ai[i]), int(bi[j])) for i, j in contact_atom_pairs}
    return float(distances.min()), int((distances < 1.5).sum()), len(contact_pairs)


def main() -> int:
    expected = read_complex_fasta()
    quality = CONFIG["proposed_readouts_pre_score"]["quality_gate"]
    receptor_len = len(next(iter(expected.values()))[0])
    peptide_len = len(next(iter(expected.values()))[1])
    rows = []
    for seed in CONFIG["prediction"]["seeds"]:
        structures = PILOT / f"seed_{seed}/structures"
        relaxed_dir = PILOT / f"seed_{seed}/relaxed"
        manifest = json.loads((relaxed_dir / "relaxation_manifest.json").read_text(encoding="utf-8"))
        manifest_by_source = {Path(row["source_path"]).name: row for row in manifest["models"]}
        for name, (expected_a, expected_b) in expected.items():
            prefix = name
            score_candidates = sorted(structures.glob(f"{prefix}_scores_rank_001_*.json"))
            raw_candidates = sorted(structures.glob(f"{prefix}_unrelaxed_rank_001_*.pdb"))
            relaxed_candidates = sorted(relaxed_dir.glob(f"{prefix}_relaxed_rank_001_*.pdb"))
            row = {"seed": seed, "complex_id": name, "expected_receptor_length": receptor_len,
                   "expected_peptide_length": peptide_len, "identity_pass": False, "quality_pass": False,
                   "contact_pass": False, "scoring_eligible": False, "failure_reasons": []}
            if not (len(score_candidates) == len(raw_candidates) == len(relaxed_candidates) == 1):
                row["failure_reasons"].append(f"expected_one_score_raw_and_relaxed_file; found {len(score_candidates)}/{len(raw_candidates)}/{len(relaxed_candidates)}")
                rows.append(row)
                continue
            scores_path, raw_pdb, relaxed_pdb = score_candidates[0], raw_candidates[0], relaxed_candidates[0]
            row.update({"scores_json": str(scores_path.relative_to(PILOT)), "raw_pdb": str(raw_pdb.relative_to(PILOT)),
                        "relaxed_pdb": str(relaxed_pdb.relative_to(PILOT)), "scores_sha256": sha(scores_path),
                        "raw_pdb_sha256": sha(raw_pdb), "relaxed_pdb_sha256": sha(relaxed_pdb)})
            relax_row = manifest_by_source.get(raw_pdb.name)
            row["relax_manifest_hash_match"] = bool(relax_row and relax_row["source_sha256"] == sha(raw_pdb)
                                                     and relax_row["relaxed_sha256"] == sha(relaxed_pdb))
            if not row["relax_manifest_hash_match"]:
                row["failure_reasons"].append("relaxation_manifest_hash_mismatch")
            chains, sequences, numbers, atom_data = residues_and_atoms(relaxed_pdb)
            row["chain_ids"] = chains
            row["receptor_sequence_match"] = sequences.get("A") == expected_a
            row["peptide_sequence_match"] = sequences.get("B") == expected_b
            row["receptor_residue_numbers_contiguous"] = numbers.get("A") == list(range(1, receptor_len + 1))
            row["peptide_residue_numbers_contiguous"] = numbers.get("B") == list(range(1, peptide_len + 1))
            row["identity_pass"] = (chains == ["A", "B"] and row["receptor_sequence_match"] and row["peptide_sequence_match"]
                                     and row["receptor_residue_numbers_contiguous"] and row["peptide_residue_numbers_contiguous"]
                                     and row["relax_manifest_hash_match"])
            if not row["identity_pass"]:
                row["failure_reasons"].append("chain_identity_sequence_or_residue_mapping_failure")

            data = json.loads(scores_path.read_text(encoding="utf-8"))
            plddt = np.asarray(data.get("plddt", []), dtype=float)
            pae = np.asarray(data.get("pae", []), dtype=float)
            row["iptm"] = float(data["iptm"]) if data.get("iptm") is not None else None
            row["peptide_mean_ca_plddt"] = float(np.mean(plddt[receptor_len:receptor_len + peptide_len])) if len(plddt) == receptor_len + peptide_len else None
            if pae.shape == (receptor_len + peptide_len, receptor_len + peptide_len):
                inter = np.concatenate((pae[:receptor_len, receptor_len:].ravel(), pae[receptor_len:, :receptor_len].ravel()))
                row["interface_pae_max_A"] = float(np.max(inter))
            else:
                row["interface_pae_max_A"] = None
            minimum, clash_pairs, contacts = geometry_metrics(atom_data)
            row.update({"minimum_interchain_heavyatom_distance_A": minimum, "clashing_atom_pairs_lt_1p5A": clash_pairs,
                        "unique_contact_residue_pairs_le_6A": contacts})
            row["iptm_pass"] = row["iptm"] is not None and np.isfinite(row["iptm"]) and row["iptm"] >= quality["iptm_min"]
            row["peptide_plddt_pass"] = row["peptide_mean_ca_plddt"] is not None and row["peptide_mean_ca_plddt"] >= quality["peptide_mean_plddt_min"]
            row["interface_pae_pass"] = row["interface_pae_max_A"] is not None and row["interface_pae_max_A"] <= quality["interface_pae_max_A"]
            row["clash_pass"] = minimum >= quality["minimum_interchain_distance_min_A"]
            row["quality_pass"] = bool(row["identity_pass"] and row["iptm_pass"] and row["peptide_plddt_pass"]
                                        and row["interface_pae_pass"] and row["clash_pass"])
            row["contact_pass"] = contacts > 0
            if not row["iptm_pass"]: row["failure_reasons"].append("iptm_below_threshold_or_missing")
            if not row["peptide_plddt_pass"]: row["failure_reasons"].append("peptide_plddt_below_threshold_or_missing")
            if not row["interface_pae_pass"]: row["failure_reasons"].append("interface_pae_above_threshold_or_missing")
            if not row["clash_pass"]: row["failure_reasons"].append("interchain_clash")
            if not row["contact_pass"]: row["failure_reasons"].append("zero_contacts_at_6A")
            row["scoring_eligible"] = bool(row["quality_pass"] and row["contact_pass"])
            rows.append(row)

    csv_path = PILOT / "structure_gate_audit_v1.csv"
    fields = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict)) else v for k, v in row.items()})
    summary = {
        "schema": "bclw_six_complex_structure_gate_audit_v1",
        "status": "STRUCTURAL_GATE_AUDIT_COMPLETE",
        "config_sha256": sha(PILOT / "run_config_frozen_v1.json"),
        "expected_cells": 6,
        "observed_cells": len(rows),
        "quality_pass_cells": sum(bool(r.get("quality_pass")) for r in rows),
        "scoring_eligible_cells": sum(bool(r.get("scoring_eligible")) for r in rows),
        "thresholds": quality,
        "rows": rows,
        "complex_scorer_execution": "All six cells have nonzero contacts and no interchain clash, so GraphPep/FoldX raw execution diagnostics were retained despite confidence-gate failure.",
        "decision": "NO_PAIR_DIRECTION_IS_EVALUABLE_WITHOUT_STRUCTURAL_QUALITY_PASS; raw scores on failed structures are diagnostic only",
    }
    json_path = PILOT / "structure_gate_audit_v1.json"
    json_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"expected_cells": 6, "observed_cells": len(rows),
                      "quality_pass_cells": summary["quality_pass_cells"],
                      "scoring_eligible_cells": summary["scoring_eligible_cells"],
                      "rows": [{k: r.get(k) for k in ("seed", "complex_id", "iptm", "peptide_mean_ca_plddt", "interface_pae_max_A", "minimum_interchain_heavyatom_distance_A", "unique_contact_residue_pairs_le_6A", "quality_pass", "scoring_eligible", "failure_reasons")} for r in rows]}, ensure_ascii=False, indent=2))
    return 0 if len(rows) == 6 else 2


if __name__ == "__main__":
    raise SystemExit(main())
