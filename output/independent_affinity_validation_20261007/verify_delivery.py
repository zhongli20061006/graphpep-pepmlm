"""Cross-check source/input/delivery integrity without inference or new scores.

Uses a separate regex extraction of the original Table 2 instead of importing
the preparation parser. This is a local secondary check, not an external
scientific review and not proof of score validity.
"""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import html
import itertools
import json
import math
from pathlib import Path
import re


TASK = Path(__file__).resolve().parent
ROOT = TASK.parent.parent
CHECKS: list[dict[str, object]] = []


def check(name: str, condition: bool, detail: object) -> None:
    CHECKS.append({"check": name, "passed": bool(condition), "detail": detail})


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(name: str) -> dict:
    return json.loads((TASK / name).read_text(encoding="utf-8"))


def load_csv(name: str) -> list[dict[str, str]]:
    with (TASK / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def plain(cell: str) -> str:
    return " ".join(html.unescape(re.sub(r"<[^>]*>", "", cell)).split())


def close(a: object, b: object) -> bool:
    return math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=1e-12)


def parse_reported(value: str) -> tuple[float, float]:
    numbers = re.findall(r"[0-9]+(?:\.[0-9]+)?", value)
    if len(numbers) != 2 or "±" not in value:
        raise ValueError(value)
    return float(numbers[0]), float(numbers[1])


def main() -> None:
    if (TASK / "verification.json").exists() or (TASK / "delivery_manifest.json").exists():
        raise FileExistsError("Refuse overwrite of an already verified delivery; use a versioned audit")
    criteria = load_json("frozen_criteria_v1.json")
    snapshot = load_json("prepared_snapshot.json")
    audit = load_json("input_audit.json")
    decision = load_json("decision.json")
    candidates = load_csv("candidate_affinity.csv")
    pairs = load_csv("mutation_pairs.csv")

    source_html = (TASK / "sources/PMC8463443.html").read_text(encoding="utf-8")
    source_tables = []
    for block in re.findall(r"<table\b[^>]*>(.*?)</table>", source_html, re.S):
        parsed_rows = []
        for row in re.findall(r"<tr\b[^>]*>(.*?)</tr>", block, re.S):
            cells = [plain(c) for c in re.findall(r"<td\b[^>]*>(.*?)</td>", row, re.S)]
            if len(cells) == 9:
                sequence = re.sub(r"\s", "", cells[0])
                if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{12}", sequence):
                    parsed_rows.append([sequence, *cells[1:]])
        if any(row[0] == "TSFAEYWNLASP" for row in parsed_rows):
            source_tables.append(parsed_rows)
    check("independent_original_Table2_extraction", len(source_tables) == 1 and len(source_tables[0]) == 23,
          {"matched_tables": len(source_tables), "expected_rows": 23})
    if not source_tables:
        raise ValueError("Original Table 2 not found")
    original = {row[0]: row for row in source_tables[0]}
    extracted = load_csv("literature_table2.csv")
    source_columns = ["sequence", "mdm2_kd_reported_nM", "mdm2_inverse_kd_ratio_reported",
                      "mdm2_ki_reported_nM", "mdm2_inverse_ki_ratio_reported",
                      "mdmx_kd_reported_nM", "mdmx_inverse_kd_ratio_reported",
                      "mdmx_ki_reported_nM", "mdmx_inverse_ki_ratio_reported"]
    check("complete_source_table_not_cherry_picked_on_scores",
          len(extracted) == 23 and all([plain(row[col]) for col in source_columns] == original[row["sequence"]]
                                      for row in extracted), "23 rows; Kd and Ki preserved in separate columns")
    check("primary_affinity_values_against_original_source",
          len(candidates) == 10 and all(
              parse_reported(original[c["sequence"]][1 if c["target"] == "MDM2" else 5]) ==
              (float(c["kd_nM"]), float(c["reported_error_nM"])) for c in candidates),
          "10 receptor-peptide observations; original SPR columns only")
    expected = {
        "MDM2_L10A": (897, 33), "MDM2_T1L": (123, 10), "MDM2_A4L": (98, 5),
        "MDM2_N8A": (229, 11), "MDM2_P12Q": (132, 2),
        "MDMX_L10A": (387, 22), "MDMX_T1L": (91, 5), "MDMX_A4L": (53, 4),
        "MDMX_N8A": (52, 4), "MDMX_S11M": (106, 3),
    }
    check("selected_values_manual_numeric_crosscheck",
          {c["peptide_id"]: (float(c["kd_nM"]), float(c["reported_error_nM"])) for c in candidates} == expected,
          "No conversion of nM to pM/uM or substitution of FP Ki")
    check("unresolved_errors_and_chemistry_preserved",
          all(not c["sd_nM"] and not c["sem_nM"] and "not_explicitly_defined" in c["error_type"]
              and c["N_terminal"].startswith("unresolved") and c["C_terminal"].startswith("unresolved")
              and c["prediction_release"] == "HOLD_data_audit" for c in candidates),
          "Missing data not invented")
    check("candidate_lengths_and_uniqueness",
          len({c["sequence"] for c in candidates}) == 6 and
          all(len(c["sequence"]) == int(c["length"]) == 12 for c in candidates) and
          Counter(c["system_id"] for c in candidates) == {"MDM2_PMI_L10A_2021": 5, "MDMX_PMI_L10A_2021": 5},
          "10 target-peptide rows / 6 unique peptides / 5 each system")

    by_id = {c["peptide_id"]: c for c in candidates}
    seen_edges = set()
    pair_valid = len(pairs) == 8
    direction_counts: Counter = Counter()
    for pair in pairs:
        p, v = by_id[pair["parent_id"]], by_id[pair["variant_id"]]
        changes = [(i+1, a, b) for i, (a, b) in enumerate(zip(p["sequence"], v["sequence"])) if a != b]
        key = (pair["system_id"], frozenset((p["peptide_id"], v["peptide_id"])))
        pair_valid &= (len(changes) == 1 and p["system_id"] == v["system_id"] == pair["system_id"] and
                       pair["parent_sequence"] == p["sequence"] and pair["variant_sequence"] == v["sequence"] and
                       key not in seen_edges and pair["count_reverse_edge_again"] == "False")
        seen_edges.add(key)
        if len(changes) == 1:
            pos, before, after = changes[0]
            pair_valid &= pair["mutation"] == f"{before}{pos}{after}" and int(pair["peptide_position_1based"]) == pos
        pk, vk = float(p["kd_nM"]), float(v["kd_nM"])
        pe, ve = float(p["reported_error_nM"]), float(v["reported_error_nM"])
        low, high = (vk-2*ve)/(pk+2*pe), (vk+2*ve)/(pk-2*pe)
        pair_valid &= (pk-2*pe > 0 and vk-2*ve > 0 and close(pair["kd_variant_over_parent"], vk/pk) and
                       close(pair["experimental_log10_gain"], math.log10(pk/vk)) and
                       close(pair["reported_error_envelope_ratio_low"], low) and
                       close(pair["reported_error_envelope_ratio_high"], high) and
                       pair["error_envelope_is_confidence_interval"] == "False")
        if pair["expected_direction"] == "stronger":
            pair_valid &= high <= 1/3 and pair["contrast_origin"] == "forward"
        else:
            pair_valid &= low >= 3 and pair["contrast_origin"] == "reverse_of_published_enhancement"
        direction_counts[(pair["system_id"], pair["expected_direction"])] += 1
    check("single_point_pairs_and_conservative_error_sensitivity", pair_valid,
          "8 distinct unordered edges; reversions declared; mean +/- 2 reported error is not a CI")
    check("within_system_primary_denominators", all(n == 2 for n in direction_counts.values()) and len(direction_counts) == 4
          and criteria["primary_metric"]["per_system_denominator"] == 4, "4 dependent contrasts per receptor; no pooling override")

    old_reference = ROOT / "output/affinity_diagnostic_mdm2_20261006/reference_affinity.csv"
    with old_reference.open(encoding="utf-8-sig", newline="") as stream:
        old = list(csv.DictReader(stream))
    development = load_csv("development_set.csv")
    check("old_ten_peptides_preserved_as_development",
          len(development) == len(old) == 10 and all(
              all(d[k] == source[k] for k in source) and d["dataset_role"] == "DEVELOPMENT_USED_FOR_METHOD_SELECTION"
              and d["eligible_as_new_primary_holdout"] == "False" for d, source in zip(development, old)),
          "Original source fields unchanged; not rebranded blind test")
    check("limited_project_exposure_claim", all(not hits for hits in audit["exact_candidate_project_exposure"].values()) and
          bool(audit["exposure_scope_limit"]) and decision["candidate_independence"]["complete_project_history_exclusion_proven"] is False,
          "Scope is the enumerated historical workspace CSVs only")

    structures = load_json("structure_source_audit.json")
    by_pdb = {s["pdb_id"]: s for s in structures["structures"]}
    uniprot = {key: "".join(line.strip() for line in (TASK / f"sources/{key}.fasta").read_text().splitlines()
                           if not line.startswith(">")) for key in ("Q00987", "O15151")}
    mdm2, mdmx = uniprot["Q00987"][24:109], uniprot["O15151"][23:108]
    mdmx_crystal = by_pdb["5UML"]["entity_sequences"]["1"]
    mismatch = [(i+24, a, b) for i, (a, b) in enumerate(zip(mdmx, mdmx_crystal)) if a != b]
    check("MDM2_sequence_and_MDMX_WT_AAA_gap", len(mdm2) == len(mdmx) == 85 and
          by_pdb["5UMM"]["entity_sequences"]["1"] == mdm2 and mismatch == [(68,"Q","A"),(69,"Q","A"),(70,"E","A")]
          and by_pdb["3EQY"]["entity_sequences"]["1"] == mdmx_crystal and audit["MDMX_exact_assay_receptor_resolved"] is False,
          {"MDM2_length": len(mdm2), "MDMX_mismatches": mismatch})
    cif_source_consistent = True
    for pdb, s in by_pdb.items():
        raw = (TASK / f"sources/{pdb}.cif").read_text()
        raw_no_space = re.sub(r"\s", "", raw)
        cif_source_consistent &= all(sequence in raw_no_space for sequence in s["entity_sequences"].values())
    check("deposited_entity_sequences_present_in_original_CIF", cif_source_consistent,
          "5 source CIFs; complete deposited sequences retained")
    check("related_crystals_not_used_as_candidate_Kd_structures",
          all(by_pdb[pdb]["entity_sequences"]["2"] == "LTFLEYWAQLMQ" for pdb in ("5UMM","5UML")) and
          all(c["exact_candidate_PDB"] == "not_identified" for c in candidates) and
          not any(c["sequence"] == "LTFLEYWAQLMQ" for c in candidates), "M3 is not a primary candidate")
    overlap = load_json("graphpep_overlap_audit.json")
    check("GraphPep_overlap_disclosure", len(overlap["pdb_id_matches"]["5UMM"]) == 2 and
          overlap["pdb_id_matches"]["5UMM"][1]["rows"] == 21 and overlap["training_clean_claim_allowed"] is False and
          decision["candidate_independence"]["not_GraphPep_training_clean"] is True,
          "5umm_B training exposure prevents a training-clean claim")

    jobs = load_csv("planned_jobs.csv")
    expected_keys = set(itertools.product(by_id, (20261021,20261022), (1,2,3)))
    actual_keys = {(job["peptide_id"],int(job["folding_seed"]),int(job["parameter_model"])) for job in jobs}
    check("budget_seeds_models_and_unstarted_jobs", len(jobs) == 60 and actual_keys == expected_keys and
          all(job["state"] == "PLANNED_NOT_STARTED_HOLD_DATA_AUDIT" for job in jobs) and
          criteria["budget"]["maximum_final_structures"] == 60,
          "Two declared seeds x three models x ten target-peptide rows; not a completed experiment")
    pdb_outputs = [str(p.relative_to(TASK)) for p in TASK.rglob("*.pdb")]
    check("no_new_predictions_or_scores", audit["actually_predicted_structures"] == audit["actually_scored_new_structures"] == 0 and
          decision["new_predictions"] == decision["new_scores"] == 0 and not pdb_outputs and
          criteria["new_scores_observed"] is False and snapshot["new_scores_observed"] is False,
          {"structure_outputs": pdb_outputs})
    check("criteria_not_modified_since_preparation", digest(TASK / "frozen_criteria_v1.json") == snapshot["criteria_sha256"],
          snapshot["criteria_sha256"])
    mismatched_files = [name for name, sha in snapshot["prepared_files"].items() if digest(TASK/name) != sha]
    mismatched_sources = [name for name, sha in snapshot["source_files"].items() if digest(TASK/"sources"/name) != sha]
    check("preparation_snapshot_files_unchanged", not mismatched_files and not mismatched_sources,
          {"mismatched_prepared": mismatched_files, "mismatched_sources": mismatched_sources})
    changed_old = [item["path"] for item in snapshot["historical_source_hashes"] if digest(ROOT/item["path"]) != item["sha256"]]
    changed_old_csvs = [item["path"] for item in audit["historical_sequence_csvs_inspected"] if digest(ROOT/item["path"]) != item["sha256"]]
    check("historical_evidence_unchanged", not changed_old and not changed_old_csvs,
          {"protected_sources": len(snapshot["historical_source_hashes"]), "CSV_sources": len(audit["historical_sequence_csvs_inspected"]),
           "changed": changed_old + changed_old_csvs})
    source_manifest = load_csv("source_manifest.csv")
    check("source_manifest_hashes", all(digest(TASK/item["file"]) == item["sha256"] and
          (TASK/item["file"]).stat().st_size == int(item["bytes"]) for item in source_manifest),
          {"source_artifacts": len(source_manifest)})
    source_audit = load_csv("candidate_source_audit.csv")
    check("two_conditional_candidates_plus_rejected_sources", len(source_audit) == 6 and
          sum(row["current_decision"].startswith("HOLD_DATA_AUDIT") for row in source_audit) == 2 and
          sum(row["current_decision"].startswith("EXCLUDE_PRIMARY") for row in source_audit) == 4,
          "No rejected Ki/IC50/modified/FEP observation filled into Kd primary panel")
    check("entry_gate_failure_not_reported_as_validation_success", decision["status"] == "HOLD_DATA_AUDIT" and
          decision["scientific_entry_gate_pass"] is False and decision["start_new_predictions_now"] is False and
          decision["not_observed"]["limited_validation_pass"] is None and decision["old_P27A_continuation_pass"] is False,
          "Data integrity and scientific entry are different outcomes")

    result = {"audit_kind": "Local secondary source/input/delivery cross-check, not external review or scoring validation",
              "checked_at_utc": datetime.now(timezone.utc).isoformat(), "checks": CHECKS,
              "delivery_integrity_pass": all(c["passed"] for c in CHECKS),
              "scientific_entry_gate_pass": False, "limited_validation_pass": None,
              "new_predictions": 0, "new_scores": 0,
              "scope_limit": "No exhaustive hidden history, homolog or all model-training audit; no experimental affinity or free-energy validation"}
    if not result["delivery_integrity_pass"]:
        failure_path = TASK / ("verification_failed_" + datetime.now(timezone.utc).strftime("%H%M%S") + ".json")
        failure_path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        print(json.dumps({"failed": [c for c in CHECKS if not c["passed"]]}, ensure_ascii=False, indent=2))
        raise SystemExit(1)
    (TASK / "verification.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    workplan = TASK / "research_workplan.md"
    workplan.write_text(workplan.read_text(encoding="utf-8").replace(
        "- [ ] 运行只读交付核验", "- [x] 运行只读交付核验"), encoding="utf-8")
    files = {p.relative_to(TASK).as_posix(): digest(p) for p in sorted(TASK.rglob("*"))
             if p.is_file() and p.name != "delivery_manifest.json" and "__pycache__" not in p.parts}
    manifest = {"state": "final_source_input_protocol_delivery_no_new_scores", "timestamp_utc": datetime.now(timezone.utc).isoformat(),
                "criteria_sha256": snapshot["criteria_sha256"], "new_predictions": 0, "new_scores": 0,
                "files_sha256": files, "scientific_entry_gate_pass": False}
    (TASK / "delivery_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"checks_passed": len(CHECKS), "delivery_integrity_pass": True,
                      "scientific_entry_gate_pass": False, "new_predictions": 0, "new_scores": 0,
                      "delivery_files": len(files)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
