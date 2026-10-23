"""Read native structures in an existing GraphPep installation; never score."""

from __future__ import annotations

import hashlib
import csv
import json
from difflib import SequenceMatcher
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
RESIDUES = dict(zip(
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(),
    "ACDEFGHIKLMNPQRSTVWY", strict=True,
))
RESIDUES["MSE"] = "M"


def sequences(path: Path) -> dict[str, str]:
    result: dict[str, list[str]] = {}
    seen = set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")) or line[17:20].strip() not in RESIDUES:
            continue
        key = (line[21:22], line[22:27])
        if key in seen:
            continue
        seen.add(key)
        result.setdefault(key[0], []).append(RESIDUES[line[17:20].strip()])
    return {chain: "".join(sequence) for chain, sequence in result.items()}


def main() -> None:
    manifest_path = PROJECT / "output/independent_direction_audit_20261020/pdz_overlap_sources/structure_input_hashes.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    native_rows = [
        row for row in manifest
        if Path(row["path"]).name in {"receptor.pdb", "peptide.pdb", "recb.pdb"}
    ]
    native_paths = [Path(row["path"]) for row in native_rows]
    peptides = {row["name"]: row["sequence"] for row in json.loads((ROOT / "peptide_sequences_source.json").read_text(encoding="utf-8"))["records"]}
    receptors = {name: row["native_core"] for name, row in json.loads((ROOT / "receptor_core_audit.json").read_text(encoding="utf-8")).items()}
    receptor_hits = []
    peptide_hits = []
    missing = []
    hash_mismatches = []
    lists = []
    base = Path("/home/g/GraphPep_v1.3")
    split_targets: dict[str, set[str]] = {"train": set(), "valid": set()}
    for split in split_targets:
        split_path = base / "Training_set" / "split_data" / f"{split}.csv"
        with split_path.open(newline="", encoding="utf-8-sig") as handle:
            split_targets[split] = {row["TARGET"].strip().lower() for row in csv.DictReader(handle)}
    list_paths = [base/"Training_set/list", base/"Training_set/split_data/train.csv", base/"Training_set/split_data/valid.csv", base/"Local_62_docking_decoys/list", base/"ADCP99_decoys/list", base/"LEADS-PEP_docking_decoys/list", base/"nr_epitope_minus_AF3_decoys/list"]
    for path in list_paths:
        payload = path.read_bytes()
        lines = payload.decode("utf-8", errors="replace").splitlines()
        matches = [line for line in lines if any(pdb in line.lower() for pdb in ("1zy3", "3fdl", "1o0l"))]
        lists.append({"path": str(path), "sha256": hashlib.sha256(payload).hexdigest(), "reference_id_matches": matches[:40], "reference_id_match_count": len(matches)})
    for row, path in zip(native_rows, native_paths, strict=True):
        if not path.is_file():
            missing.append(str(path))
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != row.get("sha256"):
            hash_mismatches.append({"path": str(path), "manifest_sha256": row.get("sha256"), "current_sha256": digest})
        group, target = path.parts[-3], path.parts[-2]
        split_membership = [split for split, targets in split_targets.items() if target.lower() in targets]
        for chain, sequence in sequences(path).items():
            if path.name == "peptide.pdb":
                for name, query in peptides.items():
                    relation = "exact" if query == sequence else "query_in_native" if query in sequence else "native_fragment_in_query" if len(sequence) >= 6 and sequence in query else None
                    if relation:
                        peptide_hits.append({"group": group, "target": target, "chain": chain, "split_membership": split_membership, "query": name, "relation": relation, "native_sequence": sequence})
            else:
                for name, query in receptors.items():
                    longest = SequenceMatcher(None, query, sequence, autojunk=False).find_longest_match(0, len(query), 0, len(sequence)).size
                    if longest >= 70 or query == sequence or query in sequence:
                        receptor_hits.append({"group": group, "target": target, "chain": chain, "split_membership": split_membership, "query": name, "native_sequence": sequence, "exact_core": query == sequence, "core_in_native": query in sequence, "native_in_core": sequence in query, "longest_exact_contiguous_block": longest})
    result = {
        "scope": "existing distributed GraphPep native structures and listed split/reference IDs only; no model inference or weights inspected",
        "native_files_examined": len(native_paths)-len(missing),
        "native_path_name_counts": {name: sum(Path(row["path"]).name == name for row in native_rows) for name in ("receptor.pdb", "recb.pdb", "peptide.pdb")},
        "missing_paths": missing,
        "manifest_sha256_mismatches": hash_mismatches,
        "split_target_counts": {split: len(targets) for split, targets in split_targets.items()},
        "source_manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "query_peptides": peptides,
        "receptor_hits_rule": "longest exact contiguous block >=70 aa, or exact/core-in-native; these are screening hits, not an exhaustive homology search",
        "receptor_hits": receptor_hits,
        "peptide_hits": peptide_hits,
        "lists": lists,
        "limits": "Core-sequence overlap is not affinity-label leakage. Lack of exact hit does not prove complete checkpoint or pretraining independence. Exact c-myc junction and reporter chemistry are not searched as full assay species because not fully represented.",
    }
    (ROOT/"graphpep_overlap_audit.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"native_files_examined":result["native_files_examined"],"missing":len(missing),"hash_mismatches":len(hash_mismatches),"recb_files":result["native_path_name_counts"]["recb.pdb"],"receptor_hits":len(receptor_hits),"peptide_hits":len(peptide_hits)},ensure_ascii=False))


if __name__ == "__main__":
    main()
