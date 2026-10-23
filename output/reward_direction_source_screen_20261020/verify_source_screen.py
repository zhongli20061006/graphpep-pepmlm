"""Verify provenance and decisive counts in the 2026-10-20 source screen."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    raw_path = ROOT / "raw" / "skempi_v2.csv"
    screen_path = ROOT / "skempi_index_screen.json"
    decision_path = ROOT / "decision.json"
    criteria_path = PROJECT / "output" / "reward_direction_benchmark_20261008" / "frozen_direction_criteria_v2.json"

    with raw_path.open("r", encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter=";"))
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    decision = json.loads(decision_path.read_text(encoding="utf-8-sig"))

    atxn_rows = [row for row in rows if row.get("#Pdb") == "4J2L_A_CD"]
    peptide_rows = [
        row for row in atxn_rows
        if row.get("Mutation(s)_PDB") in {"WC37A", "LC40S"}
    ]
    ratios = {
        row["Mutation(s)_PDB"]: float(row["Affinity_mut (M)"]) / float(row["Affinity_wt (M)"])
        for row in peptide_rows
    }
    criterion_hash = sha256(criteria_path)

    checks = {
        "raw_csv_rows_7085": len(rows) == 7085,
        "screen_groups_514": screen.get("unique_group_count") == 514,
        "screen_both_point_direction_groups_125": screen.get("groups_with_both_point_directions") == 125,
        "decision_holds_no_pilot": decision.get("status") == "HOLD_NO_PILOT",
        "frozen_criteria_hash_matches": criterion_hash == decision["frozen_criteria"]["sha256"],
        "atxn1_rows_nine": len(atxn_rows) == 9,
        "cic21_peptide_mutants_two": len(peptide_rows) == 2,
        "cic21_peptide_mutants_both_losses": len(ratios) == 2 and all(value >= 2 for value in ratios.values()),
        "new_prediction_and_scoring_counts_zero": all(value == 0 for value in decision["this_turn_compute"].values()),
    }

    artifacts = [
        ROOT / "report_zh.md",
        ROOT / "candidate_source_ranking.csv",
        decision_path,
        screen_path,
        ROOT / "screen_skempi_index.py",
        ROOT / "verify_source_screen.py",
        raw_path,
        criteria_path,
    ]
    result = {
        "as_of": "2026-10-20",
        "all_checks_pass": all(checks.values()),
        "checks": checks,
        "raw_rows": len(rows),
        "screen_groups": screen.get("unique_group_count"),
        "groups_with_both_point_directions": screen.get("groups_with_both_point_directions"),
        "atxn1_skempi_group_rows": len(atxn_rows),
        "atxn1_cic21_peptide_variant_to_wt_point_ratios": ratios,
        "frozen_criteria_sha256": criterion_hash,
        "artifacts_sha256": {
            str(path.relative_to(PROJECT)).replace("\\", "/"): sha256(path)
            for path in artifacts
        },
        "new_structure_predictions": 0,
        "new_scorer_runs": 0,
        "new_model_training_runs": 0,
        "note": "Source-screen verification only; it does not validate scorer efficacy or experimental affinity beyond the cited source records.",
    }
    output_path = ROOT / "verification.json"
    output_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"all_checks_pass": result["all_checks_pass"], "checks": checks, "ratios": ratios}, ensure_ascii=False, indent=2))
    if not result["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
