"""Independent delivery checks from CSV/A3M; does not import msa_core or models."""
import csv
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def csv_rows(name: str) -> list[dict]:
    with (OUT / name).open(encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    checks = {}
    all_rows, selected = csv_rows("candidate_decisions.csv"), csv_rows("selected_homologs.csv")
    with (ROOT / "output/bclw_construct_msa_audit_20261022/sources/uniprot_bcl2l2_candidates.tsv").open() as stream:
        source = {row["Entry"]: row for row in csv.DictReader(stream, delimiter="\t")}
    reference = load(ROOT / "output/structure_input_audit_bclw_20261020/sources/Q92843_current.json")["sequence"]["value"]
    query = reference[:164]
    checks["all_93_accessions_accounted_once"] = len(all_rows) == len(source) == 93 and {r["accession"] for r in all_rows} == set(source)
    checks["selected_csv_matches_decision_labels"] = {r["accession"] for r in selected} == {r["accession"] for r in all_rows if r["decision"] == "SELECTED"}
    checks["one_per_taxon_and_unique_core"] = len({r["taxon"] for r in selected}) == len({r["match"] for r in selected}) == len(selected)
    checks["query_identical_cores_excluded"] = all(row["match"] != query for row in selected)
    checks["full_sequences_identical_to_original_source"] = all(row["full_sequence"] == source[row["accession"]]["Sequence"] for row in all_rows)
    alignment = {r["accession"]: r for r in load(OUT / "full_alignment_evidence.json")}
    mapping_ok = True
    for row in selected:
        positions = alignment[row["accession"]]["source_positions"]
        rebuilt = "".join("-" if position is None else source[row["accession"]]["Sequence"][position - 1] for position in positions)
        mapping_ok &= len(positions) == 164 and rebuilt == row["match"]
        mapping_ok &= "".join(c for c in row["a3m"] if not c.islower()) == rebuilt
        mapping_ok &= len(row["match"].replace("-", "")) / 164 == float(row["coverage"])
        mapping_ok &= sum(a == b and a != "-" for a, b in zip(query, rebuilt)) / 164 == float(row["identity_to_query"])
    checks["source_positions_match_projection_identity_and_coverage"] = bool(mapping_ok)
    checks["selected_rows_pass_inherited_coverage_identity"] = all(float(r["coverage"]) >= 0.9 and float(r["identity_to_query"]) >= 0.4 for r in selected)
    rows = [query] + [row["match"] for row in selected]
    neighbors = [[sum(a == b and a != "-" for a, b in zip(first, second)) >= 0.8 * 164 for second in rows] for first in rows]
    exact_neff = sum(Fraction(1, sum(n)) for n in neighbors)
    representatives = []
    for i in range(len(rows)):
        if not any(neighbors[i][previous] for previous in representatives):
            representatives.append(i)
    occupied = [sum(row[column] != "-" for row in rows) / len(rows) for column in range(164)]
    values = load(OUT / "diversity_and_readiness.json")
    primary = values["primary_query_inclusive"]
    checks["independent_fraction_neff_agrees"] = math.isclose(float(exact_neff), primary["neff_80"], abs_tol=1e-12)
    checks["independent_clusters_agree"] = len(representatives) == primary["greedy_clusters_80"] and representatives == primary["greedy_representative_indices"]
    checks["all_column_occupancies_agree"] = occupied == primary["column_nongap_fraction"]
    independent_gates = {"nonquery_rows": len(selected) >= 8, "greedy_clusters": len(representatives) >= 3,
        "neff": exact_neff >= 3, "column_nongap": min(occupied) >= 0.8}
    checks["independent_frozen_gate_conjunction_agrees"] = independent_gates == values["gates"]
    decision = load(OUT / "decision.json")
    checks["decision_counts_agree"] = decision["candidate_decision_counts"] == dict(Counter(r["decision"] for r in all_rows))
    a3m = (OUT / "prepared_inputs/BCLW_Bim_A2eT_parent.a3m").read_text().splitlines()
    sequence_lines = a3m[2::2]
    checks["complex_a3m_fixed_width_and_query"] = a3m[0] == "#164,23\t1,1" and sequence_lines[0] == query + "RPEIWITQELRRIGDEFNAYYAR" and all(sum(not c.islower() for c in line) == 187 for line in sequence_lines)
    checks["homologs_receptor_only_not_paired"] = sequence_lines[1] == query + "-" * 23 and all(line.endswith("-" * 23) for line in sequence_lines[2:-1]) and sequence_lines[-1] == "-" * 164 + "RPEIWITQELRRIGDEFNAYYAR"
    checks["installed_parser_and_io_all_pass"] = load(OUT / "parser_result.json")["all_pass"] and load(OUT / "installed_io_validation.json")["all_pass"]
    sensitivity = load(OUT / "alignment_sensitivity.json")
    checks["all_selected_optima_checked_without_one_row_gate_loss"] = sensitivity["all_available_selected_optima_enumerated"] and not sensitivity["one_row_replacement_floor_failures"]
    original = load(ROOT / "output/bclw_parent_kernel_ablation_v1/prepared_protocol.json")
    protocol = load(OUT / "prepared_diagnostic_protocol.json")
    checks["quality_and_auxiliary_gates_unchanged"] = protocol["quality_gate_unchanged"] == original["quality_gate_unchanged"] and protocol["additional_diagnostic_gate"] == original["additional_diagnostic_gate"]
    command_ok = len(protocol["commands"]) == 2
    for prepared, old in zip(protocol["commands"], original["commands"]):
        command_ok &= prepared["seed"] == old["seed"]
        command_ok &= len(prepared["argv"]) == len(old["argv"])
        mode_index = old["argv"].index("--msa-mode") + 1
        differences = {i for i, (a, b) in enumerate(zip(prepared["argv"], old["argv"])) if a != b}
        command_ok &= differences == {1, 2, mode_index}
        command_ok &= prepared["argv"][mode_index] == "mmseqs2_uniref_env"
        command_ok &= "--no-use-fast-kernels" in prepared["argv"] and "--templates" not in prepared["argv"]
        command_ok &= prepared["relax_argv"][:-2] == old["relax_argv"][:-2]
    checks["future_commands_only_input_output_msa_mode_changed"] = bool(command_ok)
    checks["future_protocol_input_hash_and_freeze_match"] = sha(Path(protocol["input_a3m"])) == protocol["input_a3m_sha256"] and sha(OUT / "prepared_diagnostic_protocol.json") == (OUT / "prepared_diagnostic_protocol.sha256").read_text().strip()
    checks["implementation_freeze_unchanged_since_real_alignment"] = sha(OUT / "selection_implementation_frozen.json") == (OUT / "selection_implementation_frozen.sha256").read_text().strip()
    old_files = load(OUT / "preservation_before.json")["files"]
    changed = [row["path"] for row in old_files if not Path(row["path"]).is_file() or sha(Path(row["path"])) != row["sha256"]]
    checks["all_335_original_files_preserved"] = len(old_files) == 335 and not changed
    checks["future_output_root_absent_no_predictions_started"] = not Path(protocol["planned_output_root"]).exists() and load(OUT / "prepared_protocol_checks.json")["commands_launched"] == 0
    checks["zero_model_score_training_download_jobs"] = all(decision[key] == 0 for key in ("new_prediction_jobs", "new_relaxation_jobs", "new_reward_calls", "new_downloads", "new_training"))
    checks["new_output_under_16MiB"] = sum(p.stat().st_size for p in OUT.rglob("*") if p.is_file()) < 16777216
    result = {"checks": checks, "passed": sum(checks.values()), "total": len(checks),
        "all_pass": all(checks.values()), "old_file_mismatches": changed,
        "independent_neff_fraction": str(exact_neff), "independent_neff": float(exact_neff),
        "input_floor_gates": independent_gates,
        "scope": "Independent calculation from source CSV/A3M and original-file hashes; same investigator, not an independent biological test or external review."}
    with (OUT / "independent_verification.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    print(json.dumps(result, ensure_ascii=False))
    if not result["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
