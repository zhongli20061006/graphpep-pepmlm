"""Independently verify FoldX preparation against the original AFM coordinates.

This script does not import the benchmark driver and does not validate an
affinity score. It can run before a licensed FoldX executable is available.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "output/affinity_diagnostic_mdm2_20261006"
OUTPUT = ROOT / "output/foldx_affinity_mdm2_20261006"
AMINO_ACIDS = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))


def require(condition: bool, message: str) -> None:
    """Fail reliably even when Python is invoked with optimization enabled."""
    if not condition:
        raise ValueError(message)


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def records(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def heavy_atoms(path: Path) -> tuple[dict, dict[str, str], int]:
    atoms: dict = {}
    residues: dict[str, dict[int, str]] = {}
    hydrogens = 0
    for line in path.read_text().splitlines():
        require(not line.startswith(("HETATM", "MODEL ", "ENDMDL")),
                f"Nonstandard chemistry/model delimiters: {path}")
        if not line.startswith("ATOM  "):
            continue
        atom = line[12:16].strip()
        element = line[76:78].strip().upper()
        if not element:
            element = atom.lstrip("0123456789")[:1].upper()
        if element in {"H", "D"}:
            hydrogens += 1
            continue
        chain, residue_id = line[21], int(line[22:26])
        require(chain in {"A", "B"}, f"Unexpected chain: {path}: {chain}")
        require(line[16] == " " and line[26] == " ", f"Unresolved alternative/insertion: {path}")
        residue_name = line[17:20].strip()
        require(residue_name in AMINO_ACIDS, f"Nonstandard residue: {path}: {residue_name}")
        coordinates = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        require(all(math.isfinite(value) for value in coordinates), f"Nonfinite coordinates: {path}")
        atom_key = (chain, residue_id, atom)
        require(atom_key not in atoms, f"Duplicate atom: {path}: {atom_key}")
        atoms[atom_key] = (residue_name, coordinates, element)
        names = residues.setdefault(chain, {})
        require(residue_id not in names or names[residue_id] == residue_name, f"Mixed residue identity: {path}")
        names[residue_id] = residue_name
    require(set(residues) == {"A", "B"}, f"Missing A/B chain: {path}")
    sequences = {}
    for chain, names in residues.items():
        positions = sorted(names)
        require(positions == list(range(1, len(positions) + 1)), f"Noncontiguous numbering: {path}: {chain}")
        for position in positions:
            require(all((chain, position, atom) in atoms for atom in ("N", "CA", "C", "O")),
                    f"Incomplete backbone: {path}: {chain}{position}")
        require((chain, positions[-1], "OXT") in atoms, f"Missing free C-terminal OXT: {path}: {chain}")
        sequences[chain] = "".join(AMINO_ACIDS[names[position]] for position in positions)
    return atoms, sequences, hydrogens


def config_values(path: Path) -> dict[str, str]:
    values = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            key, value = stripped.split("=", 1)
            require(key.strip() not in values, f"Duplicate config option: {path}")
            values[key.strip()] = value.strip()
    return values


def main() -> None:
    references = records(SOURCE / "reference_affinity.csv")
    reference = {row["peptide_id"]: row for row in references}
    source_protocol = json.loads((SOURCE / "protocol.json").read_text())
    sources = records(SOURCE / "analysis/model_results.csv")
    manifest = records(OUTPUT / "input_manifest.csv")
    require(len(reference) == len(references) == 10, "Need ten unique experimental references")
    require(len(sources) == len(manifest) == 60, "Need exactly sixty source/prepared models")
    identity = lambda row: (row["peptide_id"], int(row["seed"]), int(row["model_number"]))
    indexed = {identity(row): row for row in sources}
    expected_identities = {(peptide, seed, model) for peptide in reference
                           for seed in (20261011, 20261012) for model in (1, 2, 3)}
    require(set(indexed) == expected_identities and len(indexed) == 60, "Unexpected source model identities")
    require({identity(row) for row in manifest} == expected_identities, "Unexpected prepared model identities")
    removed_hydrogens = 0
    for row in manifest:
        source = indexed[identity(row)]
        peptide = reference[row["peptide_id"]]
        source_path, input_path = Path(row["source_path"]), Path(row["input_path"])
        require(source_path == Path(source["complex_path"]), "Manifest must point to the original source")
        require(digest(source_path) == row["source_sha256"] == source["structure_sha256"],
                f"Original structure changed: {source_path}")
        require(digest(input_path) == row["input_sha256"], f"Prepared structure changed: {input_path}")
        original_atoms, original_sequences, original_hydrogens = heavy_atoms(source_path)
        prepared_atoms, prepared_sequences, prepared_hydrogens = heavy_atoms(input_path)
        require(original_atoms == prepared_atoms, f"Heavy atoms or coordinates changed: {input_path}")
        require(original_sequences == prepared_sequences, f"Sequence changed: {input_path}")
        require(prepared_sequences == {"A": source_protocol["receptor_sequence"], "B": peptide["sequence"]},
                f"Unexpected receptor/peptide sequence: {input_path}")
        require(len(prepared_sequences["A"]) == 85 and len(prepared_sequences["B"]) == 12,
                f"Incorrect chain length: {input_path}")
        require(prepared_hydrogens == 0, f"Input still contains hydrogen: {input_path}")
        require(int(row["heavy_atoms"]) == len(prepared_atoms), "Heavy atom count mismatch")
        require(int(row["removed_hydrogens"]) == original_hydrogens, "Removed hydrogen count mismatch")
        require(float(row["kd_nM"]) == float(peptide["kd_nM"]), "Experimental Kd changed")
        require(row["sequence"] == peptide["sequence"] and row["parent_id"] == peptide["parent_id"],
                "Experimental identity changed")
        quality = source_protocol["quality_gate"]
        metrics = [float(source[column]) for column in
                   ("iptm", "peptide_plddt", "interface_pae", "minimum_heavyatom_distance")]
        passing = (all(math.isfinite(value) for value in metrics)
                   and metrics[0] >= quality["iptm_min"]
                   and metrics[1] >= quality["peptide_plddt_min"]
                   and metrics[2] <= quality["interface_pae_max"]
                   and metrics[3] >= quality["minimum_distance_min"]
                   and int(source["contact_pairs_6A"]) > 0)
        require((row["quality_pass"].lower() == "true") == passing, "Quality gate provenance mismatch")
        for config_name, command in (("repair", "RepairPDB"), ("analyse", "AnalyseComplex")):
            path = input_path.parent / f"{config_name}.cfg"
            require(digest(path) == row[f"{config_name}_config_sha256"], f"Config changed: {path}")
            values = config_values(path)
            require(values["command"] == command, "Wrong FoldX command")
            require(float(values["temperature"]) == 298 and float(values["pH"]) == 7.4
                    and float(values["ionStrength"]) == 0.15, "Changed calculation conditions")
            require(values["water"] == "-IGNORE" and values["screen"].lower() == "false",
                    "Changed water/console configuration")
            require(values["noNterm"].lower() == "false" and values["noCterm"].lower() == "false",
                    "Free termini treatment changed")
            if config_name == "analyse":
                require(values["analyseComplexChains"] == "A,B", "Wrong interacting chains")
                require(values["pdb"] == "input_Repair.pdb", "Must score the repaired complex")
            else:
                require(values["pdb"] == "input.pdb", "Wrong repair input")
        removed_hydrogens += original_hydrogens
    result = {
        "status": "passed",
        "scope": "prepared_inputs_only_not_affinity_validation",
        "reference_peptides": 10,
        "prepared_models": 60,
        "all_source_structure_hashes_match": True,
        "all_heavy_atom_coordinates_unchanged": True,
        "all_receptor_peptide_sequences_complete_85_plus_12": True,
        "all_configurations_match_fixed_conditions": True,
        "removed_hydrogens_total": removed_hydrogens,
        "verification_script_sha256": digest(Path(__file__)),
    }
    (OUTPUT / "independent_input_verification.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
