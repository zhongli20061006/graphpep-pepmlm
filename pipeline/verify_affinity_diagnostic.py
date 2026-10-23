"""Independent artifact checks without importing the experiment driver."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median
from tempfile import TemporaryDirectory
import shutil

from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/affinity_diagnostic_mdm2_20261006"


def records(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def close(left: float, right: float) -> bool:
    return (math.isnan(left) and math.isnan(right)) or math.isclose(left, right, abs_tol=1e-10)


def native_sequence(path: Path, chain: str) -> str:
    names = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
    residues, atoms = {}, set()
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        residue_id = int(line[22:26])
        key = (residue_id, line[12:16].strip())
        assert key not in atoms, f"Duplicate atom {path} {key}"
        atoms.add(key)
        residues[residue_id] = names[line[17:20].strip()]
    assert list(residues) == list(range(1, len(residues) + 1))
    return "".join(residues.values())


def main() -> None:
    protocol = json.loads((OUT / "protocol.json").read_text())
    reference = records(OUT / "reference_affinity.csv")
    assert len(reference) == len({row["sequence"] for row in reference}) == 10
    for row in reference:
        assert len(row["sequence"]) == 12
        assert close(float(row["kd_M"]) * 1e9, float(row["kd_nM"]))
        assert close(float(row["sd_M"]) * 1e9, float(row["sd_nM"]))
    ref = {row["peptide_id"]: row for row in reference}
    models = records(OUT / "analysis/model_results.csv")
    assert len(models) == 60
    key = lambda row: (row["peptide_id"], row["seed"], row["model_number"])
    assert len({key(row) for row in models}) == 60
    for row in models:
        with Path(row["complex_path"]).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == row["structure_sha256"]
        assert native_sequence(Path(row["complex_path"]), "A") == protocol["receptor_sequence"]
        assert native_sequence(Path(row["complex_path"]), "B") == ref[row["peptide_id"]]["sequence"]
        passing = row["structure_pass"] == "True"
        for column in ("iptm", "peptide_plddt", "interface_pae", "minimum_heavyatom_distance"):
            assert not passing or math.isfinite(float(row[column]))
        if row["graphpep_status"] == "no_contact":
            assert not passing
            continue
        confidence = [float(pair["CONFIDENT_SCORES"]) for pair in records(Path(row["contact_feedback_path"]))]
        assert len(confidence) == int(row["n_contacts"]) == int(row["contact_pairs_6A"])
        lo = sum(weight * sum(value - .00005 > cut for value in confidence) for weight, cut in ((1, .5), (2, .7), (4, .9))) / 7
        hi = sum(weight * sum(value + .00005 > cut for value in confidence) for weight, cut in ((1, .5), (2, .7), (4, .9))) / 7
        assert -math.log1p(hi) - .000051 <= float(row["binding_score"]) <= -math.log1p(lo) + .000051
    summaries = records(OUT / "analysis/sequence_summary.csv")
    assert len(summaries) == 20
    for row in summaries:
        group = [model for model in models if (model["peptide_id"], model["seed"]) == (row["peptide_id"], row["seed"])]
        assert len(group) == 3
        passing = [model for model in group if model["structure_pass"] == "True" and model["graphpep_status"] == "ok"]
        assert int(row["passing_models"]) == len(passing)
        finite = [float(model["binding_score"]) for model in group if model["graphpep_status"] == "ok"]
        assert close(float(row["median_score_all"]), median(finite) if finite else math.nan)
        assert close(float(row["median_score_passing"]), median(float(model["binding_score"]) for model in passing) if len(passing) >= 2 else math.nan)
    indexed = {key(row): row for row in models}
    directions = records(OUT / "analysis/mutation_directions.csv")
    for row in directions:
        pairs = [(indexed[row["parent_id"], row["seed"], str(i)], indexed[row["peptide_id"], row["seed"], str(i)]) for i in (1, 2, 3)]
        gains = [float(parent["binding_score"]) - float(child["binding_score"]) for parent, child in pairs if parent["structure_pass"] == child["structure_pass"] == "True" and parent["graphpep_status"] == child["graphpep_status"] == "ok"]
        value = median(gains) if len(gains) >= 2 else math.nan
        assert int(row["common_passing_models"]) == len(gains)
        assert close(float(row["paired_median_gain"]), value)
        expected_sign = 1 if float(ref[row["peptide_id"]]["kd_nM"]) < float(ref[row["parent_id"]]["kd_nM"]) else -1
        assert (row["correct_direction"] == "True") == (math.isfinite(value) and value * expected_sign > 0)
    for row in records(OUT / "analysis/rank_correlations.csv"):
        group = [item for item in summaries if item["seed"] == row["seed"] and (row["family"] == "all" or item["family"] == row["family"]) and item[row["score_summary"]]]
        assert len(group) == int(row["n_comparable_peptides"])
        value = float(spearmanr([math.log10(float(item["kd_nM"])) for item in group], [float(item[row["score_summary"]]) for item in group]).statistic) if len(group) >= 3 else math.nan
        assert close(float(row["spearman_vs_log10_kd"]), value)
    native = records(OUT / "analysis/native_results.csv")
    common = records(OUT / "analysis/predicted_common_results.csv")
    assert len(native) == 18 and len(common) == 12
    assert sum(row["scope"] == "common" and row["pdb_id"] == "3LNZ" for row in native) == 8
    identical_inputs = {}
    for row in native:
        identical_inputs.setdefault(row["structure_sha256"], []).append(row)
    duplicate_groups = [group for group in identical_inputs.values() if len(group) > 1]
    for group in duplicate_groups:
        assert len({float(row["binding_score"]) for row in group}) == 1
    for row in [item for item in native if item["scope"] == "common"] + common:
        assert native_sequence(Path(row["complex_path"]), "A") == protocol["receptor_sequence"][2:84]
        assert native_sequence(Path(row["complex_path"]), "B") == ref[row["peptide_id"]]["sequence"][:11]
        assert int(row["contact_pairs_6A"]) == int(row["n_geometry_contact_pairs"]) == int(row["n_contacts"])
        assert len(records(Path(row["contact_feedback_path"]))) == int(row["n_contacts"])
        with Path(row["complex_path"]).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == row["structure_sha256"]
    for raw, expected in protocol["source_hashes"].items():
        with Path(raw).open("rb") as handle:
            assert hashlib.file_digest(handle, "sha256").hexdigest() == expected, raw
    # Exercise the finite-confidence gate using a copy of a real model. Original
    # coordinates and prediction confidence files remain untouched.
    from affinity_diagnostic import quality_metrics
    sample = models[0]
    with TemporaryDirectory(prefix="confidence_guard_", dir=OUT / "analysis") as temporary:
        copied_path = Path(temporary) / Path(sample["complex_path"]).name
        shutil.copyfile(sample["complex_path"], copied_path)
        values = json.loads(Path(sample["confidence_path"]).read_text())
        values["iptm"] = math.nan
        values["plddt"] = [math.nan for _ in values["plddt"]]
        copied_confidence = copied_path.with_name(copied_path.name.replace("_relaxed_", "_scores_")).with_suffix(".json")
        copied_confidence.write_text(json.dumps(values))
        metric = quality_metrics(ref[sample["peptide_id"]], int(sample["seed"]), copied_path)
        assert not metric["structure_pass"]
        assert "nonfinite_iptm" in metric["quality_fail_reasons"]
        assert "nonfinite_peptide_plddt" in metric["quality_fail_reasons"]
    result = {"status": "passed", "reference_peptides": 10, "afm_models": 60,
              "sequence_summaries_recomputed": 20, "mutation_directions_recomputed": len(directions),
              "all_crystal_copies_retained": True, "cropped_score_contact_paths_consistent": True,
              "units_and_score_direction_verified": True, "frozen_sources_verified": True,
              "nonfinite_confidence_rejected": True, "scored_structure_hashes_verified": True,
              "identical_native_input_groups_with_identical_scores": len(duplicate_groups)}
    (OUT / "analysis/independent_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
