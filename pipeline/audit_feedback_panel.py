"""Independently reconstruct the prespecified comparison from raw score CSVs."""
import argparse
import csv
from decimal import Decimal
import json
import math
from pathlib import Path
import re
from statistics import median

OUT = Path(__file__).resolve().parents[1] / "output/feedback_panel_mdm2_20261005"


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final", action="store_true")
    args = parser.parse_args()
    quality = {row["model_path"]: row["structure_pass"].lower() == "true"
               for row in read(OUT / "analysis/model_metrics.csv")}
    parent = {}
    for row in read(OUT / "parent_model_scores.csv"):
        number = int(re.search(r"_model_(\d+)_seed_", row["complex_path"]).group(1))
        parent[(row["peptide_id"], number)] = row
    children = read(OUT / "analysis/offspring_graphpep/scores.csv")
    groups = {}
    for row in children:
        groups.setdefault(row["peptide_id"], []).append(row)
    results = []
    for candidate in read(OUT / "candidates.csv"):
        if candidate["arm"] == "parent":
            continue
        group = groups.get(candidate["peptide_id"], [])
        if len(group) != 3:
            continue
        child_values, parent_values, improvements = [], [], 0
        for row in group:
            number = int(re.search(r"_model_(\d+)_seed_", row["complex_path"]).group(1))
            if number != int(row["model_number"]):
                raise ValueError("Stored model number disagrees with AFM parameter filename")
            reference = parent[(candidate["parent_id"], number)]
            if quality[row["complex_path"]] and quality[reference["complex_path"]]:
                child_score, parent_score = Decimal(row["binding_score"]), Decimal(reference["binding_score"])
                child_values.append(child_score)
                parent_values.append(parent_score)
                improvements += int(parent_score > child_score)
        gain = median(parent_values) - median(child_values) if len(child_values) >= 2 else None
        success = len(child_values) >= 2 and gain > Decimal("0.05") and improvements >= 2
        results.append({"peptide_id": candidate["peptide_id"], "parent_id": candidate["parent_id"],
                        "arm": candidate["arm"], "common_passing_models": len(child_values),
                        "median_gain_decimal": str(gain) if gain is not None else "",
                        "matched_models_improved": improvements, "predefined_improvement": bool(success)})
    if args.final:
        reported = {row["peptide_id"]: row for row in read(OUT / "analysis/sequence_summary.csv")}
        if len(results) != 36 or len(children) != 108 or len(reported) != 36:
            raise ValueError("Incomplete panel for final independent audit")
        for row in results:
            expected = reported[row["peptide_id"]]
            if row["predefined_improvement"] != (expected["predefined_improvement"].lower() == "true"):
                raise ValueError(f"Improvement flag mismatch: {row['peptide_id']}")
            if row["median_gain_decimal"]:
                if not math.isclose(float(row["median_gain_decimal"]), float(expected["median_gain_vs_predicted_parent"]), abs_tol=1e-10):
                    raise ValueError("Reported gain differs from independent decimal calculation")
            elif expected["median_gain_vs_predicted_parent"]:
                raise ValueError("Unqualified model comparison was reported as a score gain")
    if results:
        with (OUT / "analysis/independent_threshold_audit.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(results[0]))
            writer.writeheader()
            writer.writerows(results)
    summary = {}
    for row in results:
        key = f"{row['parent_id']} / {row['arm']}"
        value = summary.setdefault(key, {"children": 0, "improved": 0})
        value["children"] += 1
        value["improved"] += row["predefined_improvement"]
    print(json.dumps({"complete_child_scores": len(children), "complete_sequences": len(results),
                      "final_audit_passed": args.final, "groups": summary}, indent=2))


if __name__ == "__main__":
    main()
