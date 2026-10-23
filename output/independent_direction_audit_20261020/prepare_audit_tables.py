"""Export auditable candidate tables; these files are not prediction job inputs."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parent


def write_csv(name: str, rows: list[dict]) -> None:
    with (BASE / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    source = json.loads((BASE / "homer1_screen/homer1_single_mutants.json").read_text(encoding="utf-8"))
    rows = []
    for item in source:
        parent, variant = item["parent_sequence"], item["variant_sequence"]
        changes = [(i + 1, a, b) for i, (a, b) in enumerate(zip(parent, variant)) if a != b]
        if len(parent) != len(variant) or len(changes) != 1:
            raise ValueError("Homer1 entry is not one standard amino-acid substitution")
        position, before, after = changes[0]
        censored = item["variant_censored"]
        ratio = (item["variant_kd_lower_bound_uM"] if censored else item["variant_kd_uM"]) / item["parent_kd_uM"]
        if abs(ratio - item["ratio_point_estimate_or_lower_bound"]) > 1e-12:
            raise ValueError("Ratio does not match source records")
        rows.append({
            "background": item["background"], "source_mutation": item["mutation_name_as_source"],
            "fragment_mutation": f"{before}{position}{after}",
            "WT_fragment_sequence": parent, "variant_fragment_sequence": variant,
            "fragment_length": len(parent), "fragment_position_1based": position,
            "WT_Kd_uM": item["parent_kd_uM"], "WT_SD_uM": item["parent_sd_uM"],
            "variant_Kd_uM": "" if censored else item["variant_kd_uM"],
            "variant_SD_uM": "" if censored else item["variant_sd_uM"],
            "variant_reported_lower_bound_uM": item["variant_kd_lower_bound_uM"] if censored else "",
            "censoring": "reported_>150; detection/fitting-bound semantics require audit" if censored else "none",
            "ratio_or_lower_bound_descriptive": ratio,
            "point_effect_screen_only": item["direction_point_estimate_2fold_only"],
            "strict_primary_label": "HOLD; pairing/raw repeats and CI estimator not established",
            "replicate_n": 3, "replicate_type": "technical", "reported_error_type": "SD as reported",
            "receptor": "Homer1 EVH1; same BLI construct; not ENAH/Homer2",
            "assayed_ligand": "immobilized biotinylated SUMO-peptide; full fusion sequence/linker/site unavailable",
            "assay": "BLI association-end signal fitted as Kd; same-tip serial concentration series",
            "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC11237882/",
            "source_table": item["source_table"], "prediction_input_released": False,
        })
    write_csv("homer1_mutation_pairs.csv", rows)
    finite = [row for row in rows if row["censoring"] == "none"]
    counter = Counter(row["point_effect_screen_only"] for row in rows)
    candidates = [row for row in rows if row["background"] == "DVL3" and row["fragment_mutation"] in ("Y25F", "P21A")]
    if len(candidates) != 2:
        raise ValueError("Expected specified same-background gain/loss candidates")
    wt = candidates[0]
    pilot = [{
        "role": "WT", "fragment_id": "DVL3_WT", "fragment_sequence": wt["WT_fragment_sequence"],
        "reported_Kd_uM": wt["WT_Kd_uM"], "reported_SD_uM": wt["WT_SD_uM"],
        "future_fold_seeds_if_released": "20261008;20261009", "input_status": "HOLD_NOT_A_JOB",
        "blocking_input": "biotinylated SUMO fusion sequence/junction/site; matched representation not supplied",
    }]
    for item in candidates:
        pilot.append({
            "role": "gain_candidate" if item["fragment_mutation"] == "Y25F" else "loss_candidate",
            "fragment_id": "DVL3_" + item["fragment_mutation"], "fragment_sequence": item["variant_fragment_sequence"],
            "reported_Kd_uM": item["variant_Kd_uM"], "reported_SD_uM": item["variant_SD_uM"],
            "future_fold_seeds_if_released": "20261008;20261009", "input_status": "HOLD_NOT_A_JOB",
            "blocking_input": "same unresolved full SUMO fusion and assay-label uncertainty; fragment is not the assayed standalone ligand",
        })
    write_csv("pilot_candidate_input_audit.csv", pilot)
    candidate_audit = [
        {
            "system": "a1-syntrophin PDZ", "source_doi": "10.1038/nbt.1489",
            "primary_table": "Supplementary Table 7; PDF page 41", "raw_single_mutant_pairs": 15,
            "WT_backgrounds": 3, "gain_2fold_point_estimate": 0, "finite_loss_2fold_point_estimate": 9,
            "censored_pairs": 1, "strict_primary_labels_released": 0,
            "receptor_construct": "UPI0000566AF6 90-164, 75-aa insert; N-terminal Trx/His6 expression; exact fusion not reconstructed",
            "assayed_peptide": "5(6)-TAMRA-NNG-core10; 13 aa; C free acid inferred from explicit Wang resin",
            "assay": "fluorescence polarization with directly labeled ligand; raw three-column Kd values",
            "error_semantics": "S7 says three experimental duplicates; main caption says three replicates; biological/preparation pairing not established",
            "structure": "2PDZ reference 86-aa receptor is different from assay75 insert; reference ligand differs from full labeled13",
            "source_vs_project_exposure": "source predictor's training/test use does not imply this project's method-development use; source-discussed Oct8, labels audited Oct20; no new scores",
            "training_overlap": "separate exact/family audit; GraphPep supplied PDZ-family train35/valid4, no proved SNTA1 receptor overlap; checkpoint completeness unknown",
            "pilot_decision": "NO_GO_UNDER_FROZEN_CRITERIA", "full_benchmark_decision": "NO_GO_NO_GAIN_CLASS",
            "decisive_gap": "No >=2-fold gain; do not lower cutoff or reorient edges to manufacture class balance",
        },
        {
            "system": "Homer1 EVH1", "source_doi": "10.1002/pro.5094",
            "primary_table": "Tables 2 and 3, Homer1 BLI columns only", "raw_single_mutant_pairs": len(rows),
            "WT_backgrounds": len({row["background"] for row in rows}),
            "gain_2fold_point_estimate": counter["gain"],
            "finite_loss_2fold_point_estimate": sum(row["point_effect_screen_only"] == "loss" for row in finite),
            "censored_pairs": len(rows) - len(finite), "strict_primary_labels_released": 0,
            "receptor_construct": "134-aa expression construct in S4; TEV cleavage reported; SNA+110-aa core (113 aa) inferred, not experimentally terminal-sequenced",
            "assayed_peptide": "20/36-aa fragments fused to immobilized biotinylated SUMO; exact full fusion/linker/biotin site not supplied",
            "assay": "BLI same Homer1 EVH1; equilibrium association signal fit; proline conformer timing caveat",
            "error_semantics": "n=3 technical replicates, mean and SD; per-replicate Kd, paired structure and independent preparations unavailable",
            "structure": "1DDV reference is not the mutation series/SUMO analyte; 5ZZ9 is Homer2, not matched BLI Homer1",
            "source_vs_project_exposure": "new source-audited candidate; raw labels visible before scores, no blind/source-independent claim inferred automatically",
            "training_overlap": "targeted local exact-reference/fragment audit separate; full SUMO analyte cannot be searched as exact sequence because it is missing",
            "pilot_decision": "HOLD_INPUT_AND_LABELS; same-background DVL3 gain/loss candidates exist",
            "full_benchmark_decision": "NO_GO; at most 9 point-effect label candidates and only 2 gains",
            "decisive_gap": "exact assay ligand representation and uncertainty handling; raw11 includes 2 below threshold and 3 censored values",
        },
    ]
    write_csv("candidate_source_audit.csv", candidate_audit)
    summary = {
        "raw_pairs": len(rows), "finite_quantitative_pairs": len(finite),
        "reported_censored_pairs": len(rows) - len(finite), "point_screen": dict(counter),
        "finite_point_screen": dict(Counter(row["point_effect_screen_only"] for row in finite)),
        "maximum_label_candidates_before_uncertainty_if_all_censored_valid": counter["gain"] + counter["loss"],
        "strict_primary_labels_released": 0, "pilot_fragment_rows_not_jobs": len(pilot),
        "all_one_substitution_verified": True, "new_predictions_or_model_scores": 0,
    }
    (BASE / "homer1_numeric_audit.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
