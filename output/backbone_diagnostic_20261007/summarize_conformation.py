#!/usr/bin/env python3
"""Summarize descriptive P27A local-backbone geometry from audited CSVs."""
from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parent
ANALYSIS = ROOT / "analysis"


def read_csv(name: str) -> list[dict[str, str]]:
    with (ANALYSIS / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def median(values: list[float]) -> float:
    if not values:
        raise ValueError("cannot summarize an empty group")
    return float(statistics.median(values))


def main() -> None:
    # alpha_like is the predeclared descriptive phi/psi window; positions 8-11
    # cover the local segment ending at the P27A substitution (peptide position 11).
    torsions = [
        row
        for row in read_csv("torsions.csv")
        if row["peptide_id"] == "P53_P27A" and row["local_8_11"] == "True"
    ]
    counts: dict[tuple[str, str, str], dict[str, int]] = {}
    observed: dict[tuple[str, str, str, str], list[bool]] = {}
    for row in torsions:
        key = (row["pair_uid"], row["condition"], row["role"], row["moment"])
        observed.setdefault(key, []).append(row["alpha_like"] == "True")
    for (pair, condition, role, moment), flags in observed.items():
        if len(flags) != 4:
            raise ValueError(f"expected four local residues for {pair}/{condition}/{role}/{moment}")
        counts.setdefault((pair, condition, role), {})[moment] = sum(flags)

    expected_pairs = 30
    if len({pair for pair, _, _ in counts}) != expected_pairs:
        raise ValueError("unexpected P27A matched-pair count")
    for row in counts.values():
        if set(row) != {"before", "after"}:
            raise ValueError("incomplete before/after torsion geometry")

    alpha: dict[str, object] = {}
    for condition in ("fixed", "flexible"):
        mutant = [v for (pair, cond, role), v in counts.items() if cond == condition and role == "mutant"]
        wt_reference = [v for (pair, cond, role), v in counts.items() if cond == condition and role == "wt_reference"]
        deltas = [v["after"] - v["before"] for v in mutant]
        wt_deltas = [v["after"] - v["before"] for v in wt_reference]
        paired_changes = []
        for pair in sorted({p for p, cond, role in counts if cond == condition}):
            mut = counts[(pair, condition, "mutant")]
            wt = counts[(pair, condition, "wt_reference")]
            paired_changes.append((mut["after"] - mut["before"]) - (wt["after"] - wt["before"]))
        alpha[condition] = {
            "n_matched_pairs": len(mutant),
            "mutant_median_alpha_like_residues_8_11_before": median([float(v["before"]) for v in mutant]),
            "mutant_median_alpha_like_residues_8_11_after": median([float(v["after"]) for v in mutant]),
            "mutant_median_change_after_minus_before": median([float(v) for v in deltas]),
            "WT_reference_median_change_after_minus_before": median([float(v) for v in wt_deltas]),
            "mutant_increased_unchanged_decreased": {
                "increased": sum(v > 0 for v in deltas),
                "unchanged": sum(v == 0 for v in deltas),
                "decreased": sum(v < 0 for v in deltas),
            },
            "median_mutation_specific_change_vs_WT_reference": median([float(v) for v in paired_changes]),
            "pairs_mutant_change_exceeds_WT_reference": sum(v > 0 for v in paired_changes),
            "pairs_equal_change": sum(v == 0 for v in paired_changes),
            "pairs_mutant_change_below_WT_reference": sum(v < 0 for v in paired_changes),
        }

    # O7--N11 is the i to i+4 distance proxy whose donor is the alanine introduced
    # by P27A. It is a geometry descriptor, not a hydrogen-bond occupancy measure.
    hbond_rows = [
        row
        for row in read_csv("hbond_proxies.csv")
        if row["peptide_id"] == "P53_P27A"
        and row["role"] == "mutant"
        and row["acceptor_residue"] == "7"
        and row["donor_residue"] == "11"
    ]
    hbond: dict[str, dict[str, dict[str, float | int]]] = {}
    by_job: dict[tuple[str, str], dict[str, dict[str, str]]] = {}
    for row in hbond_rows:
        by_job.setdefault((row["job_id"], row["condition"]), {})[row["moment"]] = row
    for condition in ("fixed", "flexible"):
        rows = [v for (_, cond), v in by_job.items() if cond == condition]
        if len(rows) != expected_pairs:
            raise ValueError("unexpected P27A O7--N11 row count")
        before = [float(v["before"]["distance_A"]) for v in rows]
        after = [float(v["after"]["distance_A"]) for v in rows]
        hbond[condition] = {
            "n_structures": len(rows),
            "median_O7_N11_distance_A_before": median(before),
            "median_O7_N11_distance_A_after": median(after),
            "median_paired_distance_change_A": median([a - b for a, b in zip(after, before)]),
            "proxy_within_3_5_A_before": sum(v["before"]["proxy_present"] == "True" for v in rows),
            "proxy_within_3_5_A_after": sum(v["after"]["proxy_present"] == "True" for v in rows),
            "eligible_alanine_donors": sum(v["after"]["eligible_donor"] == "True" for v in rows),
        }

    displacement_rows = [
        row
        for row in read_csv("residue_displacements.csv")
        if row["peptide_id"] == "P53_P27A"
        and row["role"] == "mutant"
        and row["condition"] == "flexible"
        and 8 <= int(row["residue"]) <= 11
    ]
    displacement_by_job: dict[str, list[float]] = {}
    for row in displacement_rows:
        displacement_by_job.setdefault(row["job_id"], []).append(float(row["backbone_rmsd_A"]))
    per_structure_max = [max(values) for values in displacement_by_job.values()]
    if len(per_structure_max) != expected_pairs:
        raise ValueError("unexpected local displacement row count")

    result = {
        "status": "complete",
        "interpretation": "descriptive_geometry_only; not DSSP, hydrogen-bond occupancy, or affinity",
        "alpha_like_definition": "phi in [-100,-30] degrees and psi in [-80,-5] degrees",
        "segment": "peptide residues 8-11; n=30 matched WT-reference/mutant AFM model-repeat pairs per condition",
        "alpha_like": alpha,
        "O7_N11_distance_proxy": hbond,
        "flexible_local_backbone_displacement_residues_8_11": {
            "n_structures": len(per_structure_max),
            "median_per_structure_max_A": median(per_structure_max),
            "maximum_per_structure_max_A": max(per_structure_max),
            "structures_with_any_residue_over_0_01_A": sum(v > 0.01 for v in per_structure_max),
        },
    }
    (ANALYSIS / "backbone_conformation_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
