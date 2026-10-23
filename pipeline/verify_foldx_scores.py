"""Independent raw-output/statistical audit of completed FoldX screening.

No benchmark-driver functions are imported. An inputs-only preparation cannot
pass this audit; failed jobs remain explicit and never acquire an energy.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median

from scipy.stats import spearmanr

from verify_foldx_inputs import heavy_atoms

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output/foldx_affinity_mdm2_20261006"
SOURCE = ROOT / "output/affinity_diagnostic_mdm2_20261006"
SEEDS = (20261011, 20261012)
PRINCIPAL = {"PMI_N8A", "PMI_Y6A", "PMI_F3A", "PMI_W7A", "P53_P27A", "P53_L22A", "P53_L26A"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def records(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def number(value: object) -> float | None:
    if value in (None, ""):
        return None
    result = float(value)
    require(math.isfinite(result), "Nonfinite published number")
    return result


def equal_number(observed: object, expected: float | None, message: str) -> None:
    value = number(observed)
    require((value is None and expected is None) or
            (value is not None and expected is not None and math.isclose(value, expected, abs_tol=1e-10, rel_tol=1e-10)),
            message)


def boolean(value: object) -> bool | None:
    if value in (None, ""):
        return None
    require(value in (True, False, "True", "False"), f"Unexpected boolean: {value}")
    return value is True or value == "True"


def read_raw_energy(path: Path) -> float:
    """Find the explicit header and exactly one correctly identified data row."""
    require(path.name.startswith(("Interaction_", "Summary_")) and path.suffix == ".fxout",
            f"Not an interaction energy table: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        cells = [[entry.strip() for entry in row] for row in csv.reader(handle, delimiter="\t")]
    header_positions = [i for i, row in enumerate(cells) if "Interaction Energy" in row]
    require(len(header_positions) == 1, f"Missing/ambiguous raw energy header: {path}")
    position = header_positions[0]
    header = cells[position]
    require(len(set(header)) == len(header) and {"Pdb", "Group1", "Group2", "Interaction Energy"}.issubset(header),
            f"Missing/duplicate raw fields: {path}")
    data = [row for row in cells[position + 1:] if row and any(row)]
    require(len(data) == 1 and len(data[0]) == len(header), f"Ambiguous/malformed raw rows: {path}")
    row = dict(zip(header, data[0]))
    require({row["Group1"], row["Group2"]} == {"A", "B"}, f"Wrong interacting chains: {path}")
    require(Path(row["Pdb"]).stem == "input_Repair", f"Wrong scored structure: {path}")
    value = number(row["Interaction Energy"])
    require(value is not None, f"Empty raw energy: {path}")
    return value


def correlation(energies: dict[str, float], references: dict[str, dict], family: str) -> tuple[float | None, int]:
    subset = sorted(peptide for peptide in energies if family == "overall" or references[peptide]["family"] == family)
    if len(subset) < 2 or len({energies[peptide] for peptide in subset}) < 2:
        return None, len(subset)
    rho = float(spearmanr([energies[peptide] for peptide in subset],
                         [math.log10(float(references[peptide]["kd_nM"])) for peptide in subset]).statistic)
    require(math.isfinite(rho), "Nonfinite independently computed correlation")
    return rho, len(subset)


def main() -> None:
    protocol = json.loads((OUTPUT / "protocol.json").read_text())
    state = json.loads((OUTPUT / "execution_state.json").read_text())
    require(state["state"] in {"complete", "complete_with_failures"},
            "Real FoldX run is not complete; this is not an affinity-validation result")
    tool = json.loads((OUTPUT / "tool_manifest.json").read_text())
    for raw, expected in protocol["source_hashes"].items():
        require(sha(Path(raw)) == expected, f"Frozen source changed: {raw}")
    for relative, expected in protocol["prepared_hashes"].items():
        require(sha(OUTPUT / relative) == expected, f"Prepared artifact changed: {relative}")
    fingerprint = tool["fingerprint"]
    require(sha(Path(fingerprint["executable"])) == fingerprint["executable_sha256"], "FoldX executable changed")
    for resource in fingerprint["resources"].values():
        require(sha(Path(resource["source"])) == resource["sha256"], "FoldX resource changed")
    for relative, expected in tool["capture_hashes"].items():
        require(sha(OUTPUT / relative) == expected, "Tool help capture changed")
    reference_rows = records(SOURCE / "reference_affinity.csv")
    references = {row["peptide_id"]: row for row in reference_rows}
    manifest = records(OUTPUT / "input_manifest.csv")
    require(len(references) == 10 and len(manifest) == len(state["jobs"]) == 60, "Incomplete real benchmark")
    model_rows = records(OUTPUT / "analysis/model_results.csv")
    published_models = {row["job_id"]: row for row in model_rows}
    require(len(published_models) == len(model_rows) == 60, "Incomplete/duplicate published model table")
    actual: dict[str, float | None] = {}
    successful = 0
    for row in manifest:
        job = OUTPUT / "jobs" / row["job_id"]
        require(sha(job / "result.json") == state["jobs"][row["job_id"]], "Changed result record")
        result = json.loads((job / "result.json").read_text())
        artifacts = {path.relative_to(job).as_posix(): sha(path) for path in job.rglob("*")
                     if path.is_file() and path.name not in {"input.pdb", "repair.cfg", "analyse.cfg", "result.json"}}
        require(artifacts == result["artifact_hashes"], f"Changed raw artifacts: {job.name}")
        published = published_models[row["job_id"]]
        require(result["status"] == published["status"] and published["units"] == result["units"] == "kcal/mol",
                "Status/units mismatch")
        require(published["input_sha256"] == row["input_sha256"] and published["source_sha256"] == row["source_sha256"],
                "Published structure provenance mismatch")
        if result["status"] == "failed":
            require("interaction_energy_kcal_mol" not in result and not published["interaction_energy_kcal_mol"],
                    "Failed job carries an energy")
            require(result["failure_reason"] == published["failure_reason"] and bool(result["failure_reason"]),
                    "Failure reason missing/changed")
            actual[row["job_id"]] = None
            continue
        require(result["status"] == "ok", "Unknown job state")
        for stage in ("repair", "analyse"):
            command = json.loads((job / f"{stage}.command.json").read_text())
            require(command["returncode"] == 0 and not command["error"] and
                    command["command"] == [fingerprint["executable"], "-f", f"{stage}.cfg"],
                    "Successful energy lacks successful FoldX command")
        _, sequences, _ = heavy_atoms(job / "input_Repair.pdb")
        require(sequences == {"A": protocol["scientific_protocol"].get("receptor_sequence", "")
                             or json.loads((SOURCE / "protocol.json").read_text())["receptor_sequence"],
                             "B": row["sequence"]}, "Repaired chemistry/sequence changed")
        candidates = sorted(job.glob("Interaction_*.fxout")) or sorted(job.glob("Summary_*.fxout"))
        require(len(candidates) == 1, "Missing/ambiguous raw energy file")
        energy = read_raw_energy(candidates[0])
        for summary in job.glob("Summary_*.fxout"):
            equal_number(read_raw_energy(summary), energy, "Summary and interaction energy disagree")
        equal_number(result["interaction_energy_kcal_mol"], energy, "Result energy differs from raw")
        equal_number(published["interaction_energy_kcal_mol"], energy, "Published energy differs from raw")
        require(result["raw_energy_paths"] == [str(candidates[0])], "Raw result path mismatch")
        actual[row["job_id"]] = energy
        successful += 1
    summaries = records(OUTPUT / "analysis/sequence_summary.csv")
    directions = records(OUTPUT / "analysis/mutation_directions.csv")
    correlations = records(OUTPUT / "analysis/rank_correlations.csv")
    require(len(summaries) == 40 and len(directions) == 32 and len(correlations) == 12, "Unexpected analysis table counts")
    summary_map = {(row["scope"], row["peptide_id"], int(row["seed"])): row for row in summaries}
    direction_map = {(row["scope"], row["peptide_id"], int(row["seed"])): row for row in directions}
    correlation_map = {(row["scope"], int(row["seed"]), row["family"]): row for row in correlations}
    require(len(summary_map) == 40 and len(direction_map) == 32 and len(correlation_map) == 12, "Duplicate analysis identities")
    decision = json.loads((OUTPUT / "analysis/decision.json").read_text())
    require(decision["scoring_complete"] is True, "Completed raw run is not marked complete in decision")
    gate_map = {int(row["seed"]): row for row in decision["seed_gates"]}
    require(set(gate_map) == set(SEEDS) and len(decision["seed_gates"]) == 2, "Wrong decision seeds")
    independent_gates = []
    for scope in ("primary_quality", "all_models"):
        minimum = 2 if scope == "primary_quality" else 1
        for seed in SEEDS:
            selected: dict[str, dict[int, float]] = {}
            energies: dict[str, float] = {}
            for peptide in references:
                rows = [row for row in manifest if row["peptide_id"] == peptide and int(row["seed"]) == seed]
                values = {int(row["model_number"]): actual[row["job_id"]] for row in rows
                          if actual[row["job_id"]] is not None and
                          (scope == "all_models" or row["quality_pass"] == "True")}
                eligible = len(values) >= minimum
                selected[peptide] = values if eligible else {}
                value = median(values.values()) if eligible else None
                if eligible:
                    energies[peptide] = value
                published = summary_map[scope, peptide, seed]
                require(boolean(published["eligible"]) == eligible and int(published["n_scored_selected"]) == len(values),
                        "Summary eligibility/count mismatch")
                equal_number(published["median_energy_kcal_mol"], value, "Summary median mismatch")
                equal_number(published["minimum_energy_kcal_mol"], min(values.values()) if values else None, "Summary min mismatch")
                equal_number(published["maximum_energy_kcal_mol"], max(values.values()) if values else None, "Summary max mismatch")
                equal_number(published["model_range_kcal_mol"], max(values.values()) - min(values.values()) if values else None,
                             "Summary model variation mismatch")
            correct: dict[str, bool | None] = {}
            for peptide, reference in references.items():
                parent = reference["parent_id"]
                if parent == peptide:
                    continue
                shared = sorted(set(selected[parent]) & set(selected[peptide]))
                gain = median(selected[parent][model] - selected[peptide][model] for model in shared) if len(shared) >= minimum else None
                expected = 1 if float(reference["kd_nM"]) < float(references[parent]["kd_nM"]) else -1
                success = gain * expected > 0 if gain is not None else None
                published = direction_map[scope, peptide, seed]
                require(int(published["n_matched_models"]) == len(shared) and
                        published["matched_models"] == ";".join(map(str, shared)), "Paired model identities mismatch")
                require(boolean(published["comparable"]) == (gain is not None) and
                        boolean(published["direction_correct"]) == success, "Direction classification mismatch")
                require(boolean(published["primary_direction_case"]) == (peptide in PRINCIPAL), "Principal direction set changed")
                equal_number(published["gain_kcal_mol"], gain, "Paired difference median mismatch")
                correct[peptide] = success
            rhos = {}
            for family in ("overall", "PMI", "p53"):
                rho, count = correlation(energies, references, family)
                published = correlation_map[scope, seed, family]
                require(int(published["n_peptides"]) == count, "Correlation sample count mismatch")
                equal_number(published["rho_energy_log10Kd"], rho, "Rank correlation mismatch")
                rhos[family] = rho
            if scope == "primary_quality":
                all_comparable = all(correct[peptide] is not None for peptide in PRINCIPAL)
                correct_count = sum(correct[peptide] is True for peptide in PRINCIPAL)
                required_correct = all(correct[peptide] is True for peptide in ("PMI_N8A", "P53_P27A"))
                pass_gate = (len(energies) == 10 and all_comparable and correct_count >= 6 and required_correct
                             and rhos["overall"] is not None and rhos["overall"] >= 0.6
                             and all(rhos[family] is not None and rhos[family] > 0 for family in ("PMI", "p53")))
                published = gate_map[seed]
                require(published["all_ten_eligible"] == (len(energies) == 10) and
                        published["n_comparable_peptides"] == len(energies), "Incomplete panel gate mismatch")
                require(published["all_seven_comparable"] == all_comparable and published["correct_of_seven"] == correct_count
                        and published["N8A_and_P27A_correct"] == required_correct and published["gate_pass"] == pass_gate,
                        "Frozen continuation decision mismatch")
                for family in rhos:
                    equal_number(published["rho"][family], rhos[family], "Decision correlation mismatch")
                independent_gates.append(pass_gate)
    require(decision["continuation_pass"] == all(independent_gates), "Continuation state mismatch")
    require(decision["counts"] == {"ok": successful, "failed": 60 - successful, "not_run": 0}, "Real scoring count mismatch")
    result = {"status": "passed", "scope": "independent_raw_output_and_statistics_audit",
              "real_scored_models": successful, "failed_models_retained": 60 - successful,
              "sequence_summaries_recomputed": 40, "mutation_directions_recomputed": 32,
              "rank_correlations_recomputed": 12, "frozen_gates_recomputed": 2,
              "affinity_screen_pass": all(independent_gates),
              "verification_script_sha256": sha(Path(__file__)),
              "limitation": "Audit correctness does not establish physical free-energy accuracy or holdout generalization"}
    (OUTPUT / "analysis/independent_score_verification.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
