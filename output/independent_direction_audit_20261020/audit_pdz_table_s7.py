"""Curate publisher Table S7 without scoring models or relaxing label criteria."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import statistics
from pathlib import Path

import pymupdf

BASE = Path(__file__).resolve().parent
SOURCE = BASE / "sources" / "pdz_supplement_publisher.pdf"
CRITERIA = BASE.parent / "reward_direction_benchmark_20261008" / "frozen_direction_criteria_v2.json"
SOURCE_URL = (
    "https://media.springernature.com/original/springer-static/esm/"
    "art%3A10.1038%2Fnbt.1489/MediaObjects/41587_2008_BFnbt1489_MOESM436_ESM.pdf"
)
# Independently transcribed from the rendered original page, including -1 sentinel.
VISUALLY_CHECKED_ROWS = [
    ("CLDTSRETDL", ("1.6", "1.4", "1.9")),
    ("CLDTSKETDL", ("3.1", "3.0", "3.1")),
    ("CLDTSRDTDL", ("47.8", "85.3", "113.5")),
    ("CLDTSRESDL", ("5.3", "5.3", "5.4")),
    ("CLDTSRETWL", ("4.8", "4.7", "5.1")),
    ("CLDTSRETDV", ("0.9", "0.9", "1.0")),
    ("SPDRDRESIV", ("1.1", "1.0", "1.1")),
    ("SPDRDMESIV", ("5.9", "6.2", "6.2")),
    ("SPDRDRRSIV", ("34.0", "37.9", "39.4")),
    ("SPDRDRETIV", ("1.0", "0.9", "1.0")),
    ("SPDRDRESDV", ("2.2", "2.3", "2.3")),
    ("SPDRDRESIL", ("4.8", "4.4", "4.8")),
    ("NLKAGRETTV", ("1.6", "1.7", "1.7")),
    ("NLKAGTETTV", ("13.1", "11.5", "12.9")),
    ("NLKAGRSTTV", ("19.3", "19.4", "20.3")),
    ("NLKAGREYTV", ("-1", "-1", "-1")),
    ("NLKAGRETQV", ("2.4", "2.3", "2.3")),
    ("NLKAGRETTL", ("2.2", "2.1", "2.3")),
]
BACKGROUNDS = ("Kv1.5", "Nav1.5", "KIF1B")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    with pymupdf.open(SOURCE) as document:
        text = document[40].get_text()
    pattern = r"([ACDEFGHIKLMNPQRSTVWY]{10})\s*(?:\([^)]*\))?\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)"
    extracted = [(match[0], tuple(match[1:])) for match in re.findall(pattern, text)]
    if extracted != VISUALLY_CHECKED_ROWS:
        raise ValueError("PDF extraction does not match independently checked source transcription")
    if len({sequence for sequence, _ in extracted}) != 18:
        raise ValueError("Expected 18 unique source ligands")
    raw_rows, ligands, pairs = [], [], []
    for index, (sequence, raw_values) in enumerate(extracted):
        background = BACKGROUNDS[index // 6]
        wt_sequence, wt_raw = extracted[(index // 6) * 6]
        mutations = [(i + 1, a, b) for i, (a, b) in enumerate(zip(wt_sequence, sequence)) if a != b]
        is_wt = index % 6 == 0
        if len(mutations) != (0 if is_wt else 1):
            raise ValueError(f"Invalid one-substitution relationship: {sequence}")
        no_binding = raw_values == ("-1", "-1", "-1")
        if not no_binding and any(float(value) <= 0 for value in raw_values):
            raise ValueError("Non-positive values must be source noninteraction sentinels only")
        values = [] if no_binding else [float(value) for value in raw_values]
        ligand_id = f"a1syn_{background}_{sequence}"
        for column, value in enumerate(raw_values, start=1):
            raw_rows.append({
                "ligand_id": ligand_id, "background": background,
                "core_sequence_as_reported": sequence,
                "assayed_amino_acid_sequence": "NNG" + sequence,
                "source_table": "Supplementary Table 7", "pdf_page_1based": 41,
                "measurement_column": column, "raw_source_value": value,
                "raw_unit": "uM", "kd_uM": "" if no_binding else value,
                "measurement_status": "noninteraction_sentinel" if no_binding else "quantitative",
                "replicate_header_verbatim": "three experimental duplicates",
                "replicate_independence_and_pairing": "not established",
                "source_url": SOURCE_URL,
            })
        row = {
            "ligand_id": ligand_id, "background": background, "is_WT": is_wt,
            "core_sequence_as_reported": sequence, "core_length": len(sequence),
            "assayed_amino_acid_sequence": "NNG" + sequence, "assayed_aa_length": 13,
            "N_terminal_chemistry": "5(6)-TAMRA attached to N-terminal NNG linker",
            "C_terminal_chemistry": "free acid inferred from explicit Wang-resin synthesis/cleavage; not a direct terminal-group statement",
            "Kd_source_col1_uM": raw_values[0], "Kd_source_col2_uM": raw_values[1],
            "Kd_source_col3_uM": raw_values[2],
            "Kd_mean_uM_descriptive": "" if no_binding else statistics.mean(values),
            "Kd_sample_SD_uM_descriptive_only": "" if no_binding else statistics.stdev(values),
            "n_reported_columns": 3, "biological_n": "unknown",
            "censoring": "noninteraction_-1; detection limit not established" if no_binding else "none",
            "assay": "fluorescence polarization; directly labeled ligand",
            "receptor": "mouse a1-syntrophin PDZ 75-aa insert; N-terminal Trx/His6; complete tag/linker sequence unresolved",
            "receptor_range_reported_S2": "UPI0000566AF6 90-164; 75 aa; confirmed by original UniParc sequence; distinct from tutorial/PDB 86-aa reference",
            "source_url": SOURCE_URL, "source_pdf_page_1based": 41,
            "model_input_release": "HOLD; TAMRA/tag/construct representation not resolved",
        }
        ligands.append(row)
        if is_wt:
            continue
        wt_mean = statistics.mean(float(value) for value in wt_raw)
        ratio = None if no_binding else statistics.mean(values) / wt_mean
        point_class = (
            "not_classifiable_no_binding_limit_unknown" if ratio is None else
            "gain_at_least_2fold_point_estimate_only" if ratio <= 0.5 else
            "loss_at_least_2fold_point_estimate_only" if ratio >= 2 else
            "below_2fold_effect_threshold"
        )
        position, ref, alt = mutations[0]
        pairs.append({
            "pair_id": f"{background}_{ref}{position}{alt}", "background": background,
            "WT_core_sequence": wt_sequence, "variant_core_sequence": sequence,
            "WT_assayed_aa_sequence": "NNG" + wt_sequence,
            "variant_assayed_aa_sequence": "NNG" + sequence,
            "core_position_1based": position, "assayed_aa_position_1based": position + 3,
            "position_relative_C_terminus": position - 10,
            "mutation": f"{ref}{position}{alt}", "mutation_count": 1,
            "WT_Kd_mean_uM_descriptive": wt_mean,
            "variant_Kd_mean_uM_descriptive": "" if no_binding else statistics.mean(values),
            "ratio_variant_over_WT_descriptive": "" if ratio is None else ratio,
            "binding_improvement_fold_descriptive": "" if ratio is None else 1 / ratio,
            "point_effect_screen": point_class,
            "strict_primary_label": "HOLD_UNCERTAINTY_OR_CENSORING; no CI reconstructed",
            "qualifies_frozen_gain_screen": ratio is not None and ratio <= 0.5,
            "source_table": "Supplementary Table 7", "source_pdf_page_1based": 41,
        })
    write_csv(BASE / "pdz_s7_raw_measurements.csv", raw_rows)
    write_csv(BASE / "pdz_ligand_input_audit.csv", ligands)
    write_csv(BASE / "pdz_mutation_pairs.csv", pairs)
    gains = [pair for pair in pairs if pair["qualifies_frozen_gain_screen"]]
    losses = [pair for pair in pairs if pair["point_effect_screen"].startswith("loss_")]
    quantitative = [pair for pair in pairs if pair["ratio_variant_over_WT_descriptive"] != ""]
    best = min(quantitative, key=lambda pair: pair["ratio_variant_over_WT_descriptive"])
    checks = {
        "source_pdf_sha256": sha256(SOURCE), "criteria_source": str(CRITERIA),
        "criteria_sha256": sha256(CRITERIA), "source_visual_check_pages": [4, 41],
        "exact_pdf_vs_visual_transcription_match": True,
        "source_ligands": len(ligands), "WT_backgrounds": 3,
        "single_mutant_pairs": len(pairs), "raw_source_values": len(raw_rows),
        "quantitative_values": sum(row["kd_uM"] != "" for row in raw_rows),
        "noninteraction_sentinel_values": sum(row["kd_uM"] == "" for row in raw_rows),
        "all_assayed_aa_sequences_length_13": all(len(row["assayed_amino_acid_sequence"]) == 13 for row in ligands),
        "all_pairs_hamming_distance_one": all(pair["mutation_count"] == 1 for pair in pairs),
        "twofold_point_estimate_gain_pairs": len(gains),
        "twofold_point_estimate_loss_candidates": len(losses),
        "strict_primary_labels_released": 0,
        "best_point_estimate_gain": best,
        "decision": "NO_GO_PDZ_PILOT_UNCHANGED_CRITERIA",
        "decisive_reason": "No mutant achieves the frozen >=2-fold enhancement even at the point-estimate screen; uncertainty and chemistry cannot repair the absent gain class",
        "new_model_predictions_or_scores": 0,
    }
    (BASE / "pdz_numeric_audit.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in checks.items() if key != "best_point_estimate_gain"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
