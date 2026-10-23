"""Build and audit a small local MSA; never invokes a prediction or scorer."""
import ast
import csv
import hashlib
import json
import platform
import signal
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import Bio
from Bio import Align
from Bio.Align import substitution_matrices

import msa_core as core

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
BATCH = Path("/home/g/miniconda3/envs/colabfold_exp0/lib/python3.11/site-packages/colabfold/batch.py")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def save(name: str, value) -> None:
    with (OUT / name).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


def write_csv(name: str, records: list[dict]) -> None:
    fields = list(dict.fromkeys(key for record in records for key in record))
    with (OUT / name).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict)) else value
                             for key, value in record.items()})


def snapshot_originals() -> list[dict]:
    paths = set()
    original_audit = ROOT / "output/bclw_construct_msa_audit_20261022"
    paths.update(path for path in original_audit.rglob("*") if path.is_file())
    for directory in ("bh3_bclw_six_complex_pilot_20261020", "bclw_parent_kernel_ablation_v1"):
        paths.update(path for path in (ROOT / "output" / directory).rglob("*") if path.is_file())
    for name in (
        "PROJECT_HANDOFF_2026-10-20.md",
        "output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json",
        "output/structure_input_audit_bclw_20261020/sources/Q92843_current.json",
        "output/structure_input_audit_bclw_20261020/prepared_inputs/parent_only.fasta",
        "output/bh3_bclw_source_audit_20261020/receptor_core_audit.json",
        "output/receptor_context_disagreement_audit_20261022/summary.json",
        "output/receptor_context_pmi_development_comparison_20261022/summary.json",
    ):
        paths.add(ROOT / name)
    paths.add(BATCH)
    return [{"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}
            for path in sorted(paths)]


def check_old_manifests() -> dict:
    results = {}
    for directory, filename in (
        ("bh3_bclw_six_complex_pilot_20261020", "run_result_manifest_v1.json"),
        ("bclw_parent_kernel_ablation_v1", "artifact_manifest.json"),
    ):
        base = ROOT / "output" / directory
        manifest = load(base / filename)
        mismatches = []
        for item in manifest["files"]:
            path = base / item["path"]
            if not path.exists() or path.stat().st_size != item["bytes"] or digest(path) != item["sha256"]:
                mismatches.append(item["path"])
        results[directory] = {"entries": len(manifest["files"]), "mismatches": mismatches}
    return results


def parse_with_installed_function(content: str, receptor: str, peptide: str,
                                  expected_a3m: list[str]) -> dict:
    """Isolate actual parser functions; mock only template feature creation."""
    source = BATCH.read_text()
    tree = ast.parse(source)
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef)
             and node.name in {"normalize_a3m", "unserialize_msa"}]
    if len(nodes) != 2:
        raise ValueError("Installed parser definitions changed")
    snippet = "from __future__ import annotations\n\n" + "\n\n".join(ast.get_source_segment(source, node) for node in nodes)
    (OUT / "installed_parser_used.py.txt").write_text(snippet, encoding="utf-8")
    namespace = {"mk_mock_template": lambda sequence: {"mock_template_for": sequence}}
    exec(compile(snippet, str(BATCH) + "::isolated_parser", "exec"), namespace)
    unpaired, paired, query_seqs, cardinalities, _ = namespace["unserialize_msa"](
        [content], [receptor, peptide])
    receptor_rows = unpaired[0].splitlines()[1::2]
    peptide_rows = unpaired[1].splitlines()[1::2]
    paired_rows = [part.splitlines()[1::2] for part in paired]
    checks = {
        "query_sequences_exact": query_seqs == [receptor, peptide],
        "cardinalities_exact": cardinalities == [1, 1],
        "receptor_rows_exact": receptor_rows == [receptor] + expected_a3m,
        "peptide_only_query": peptide_rows == [peptide],
        "paired_rows_only_concatenated_query": paired_rows == [[receptor], [peptide]],
        "receptor_match_width_164": all(sum(not aa.islower() for aa in row) == 164 for row in receptor_rows),
        "no_peptide_insertion_leakage": all(row == peptide for row in peptide_rows),
    }
    save("parser_result.json", {
        "installed_batch_path": str(BATCH), "installed_batch_sha256": digest(BATCH),
        "scope": "Exact installed normalize_a3m/unserialize_msa function bodies, template factory stubbed; no ColabFold/JAX import or complete feature/prediction validation.",
        "checks": checks, "all_pass": all(checks.values()),
        "unpaired_msa": unpaired, "paired_msa": paired,
    })
    return checks


def main() -> None:
    started = time.monotonic()
    started_utc = datetime.now(timezone.utc).isoformat()
    config = load(OUT / "selection_implementation_frozen.json")
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError("CPU audit budget exceeded")))
    signal.alarm(config["cpu_worker_wall_seconds_max"])
    old_snapshot = snapshot_originals()
    save("preservation_before.json", {"started_utc": started_utc, "files": old_snapshot})
    (OUT / "selection_implementation_frozen.sha256").write_text(digest(OUT / "selection_implementation_frozen.json") + "\n")
    inherited_path = ROOT / config["inherited_rules_path"]
    inherited = load(inherited_path)["msa_selection_frozen_before_retrieval"]
    old_manifests = check_old_manifests()
    assert all(not value["mismatches"] for value in old_manifests.values()), old_manifests
    assert digest(ROOT / config["source_tsv_path"]) == config["source_tsv_expected_sha256"]
    assert digest(ROOT / config["parent_fasta_path"]) == config["parent_fasta_expected_sha256"]
    assert digest(ROOT / "output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json") == "23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b"
    reference_record = load(ROOT / config["reference_path"])
    assert reference_record["primaryAccession"] == "Q92843"
    reference = reference_record["sequence"]["value"]
    receptor = reference[:164]
    assert len(reference) == 193 and receptor[132] == "Q"
    original_core = load(ROOT / "output/bh3_bclw_source_audit_20261020/receptor_core_audit.json")["Bcl-w"]["native_core"]
    fasta = (ROOT / config["parent_fasta_path"]).read_text().splitlines()
    actual_receptor, peptide = "".join(fasta[1:]).split(":")
    assert receptor == original_core == actual_receptor
    assert peptide == "RPEIWITQELRRIGDEFNAYYAR" and len(peptide) == 23
    save("input_correspondence.json", {
        "reference_accession": "Q92843", "reference_length": len(reference),
        "receptor_sequence": receptor, "peptide_sequence": peptide,
        "receptor_length": len(receptor), "peptide_length": len(peptide),
        "Q133_unchanged": receptor[132] == "Q", "matches_old_fasta_and_core_audit": True,
        "model_representation": "Unchanged 1-164 core + standard parent 23mer; not exact assay chemistry.",
        "omissions": ["N-terminal fluorescein", "C-terminal amidation", "N-terminal receptor c-myc/TEV junction of unrecovered complete sequence"],
        "uncertainty": "Assay clone Q133 identity not independently recovered; do not silently change to R133.",
    })
    with (ROOT / config["source_tsv_path"]).open(encoding="utf-8-sig") as stream:
        raw = list(csv.DictReader(stream, delimiter="\t"))
    assert len(raw) == 93
    assert len({row["Entry"] for row in raw}) == 93
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.epsilon = 1e-6
    (OUT / "aligner_parameters.txt").write_text(str(aligner), encoding="utf-8")
    records, alignment_evidence = [], []
    for row in raw:
        record = {"accession": row["Entry"], "entry_name": row["Entry Name"],
                  "reviewed": row["Reviewed"].lower() == "reviewed", "taxon": row["Organism (ID)"],
                  "organism": row["Organism"], "protein_name": row["Protein names"],
                  "gene_names": row["Gene Names"], "full_length": len(row["Sequence"]),
                  "full_sequence": row["Sequence"], "source_annotation_score": row["Annotation"]}
        reason = core.metadata_reason(row, *inherited["full_length_inclusive"])
        if reason:
            record["decision"] = reason
            records.append(record)
            continue
        alignments = aligner.align(reference, row["Sequence"])
        alignment = alignments[0]
        try:
            optimal_count = len(alignments)
        except OverflowError:
            optimal_count = "OVERFLOW_MANY_OPTIMA"
        projection = core.project_core(reference, row["Sequence"], alignment.coordinates.tolist(), 164)
        coverage = sum(aa != "-" for aa in projection["match"]) / 164
        identity = core.identity(receptor, projection["match"])
        record.update({"coverage": coverage, "identity_to_query": identity,
                       "match": projection["match"], "a3m": projection["a3m"],
                       "alignment_score": alignment.score, "optimal_alignment_count": optimal_count})
        alignment_evidence.append({"accession": row["Entry"], "full_alignment": str(alignment),
                                   "coordinates": alignment.coordinates.tolist(), **projection})
        if coverage < inherited["core_coverage_min"]:
            record["decision"] = "EXCLUDE_CORE_COVERAGE"
        elif identity < inherited["core_matches_divided_by_164_min"]:
            record["decision"] = "EXCLUDE_CORE_IDENTITY"
        elif projection["match"] == receptor:
            record["decision"] = "EXCLUDE_QUERY_IDENTICAL_CORE"
        else:
            record["decision"] = "ELIGIBLE_BEFORE_REPRESENTATIVE_SELECTION"
        records.append(record)
    seen_taxa, seen_cores = {}, {}
    finalists = []
    eligible = sorted([r for r in records if r["decision"] == "ELIGIBLE_BEFORE_REPRESENTATIVE_SELECTION"], key=core.priority)
    for record in eligible:
        if record["taxon"] in seen_taxa:
            record.update(decision="EXCLUDE_SAME_TAXON", representative=seen_taxa[record["taxon"]])
        else:
            seen_taxa[record["taxon"]] = record["accession"]
            finalists.append(record)
    unique = []
    for record in finalists:
        if record["match"] in seen_cores:
            record.update(decision="EXCLUDE_DUPLICATE_CORE", representative=seen_cores[record["match"]])
        else:
            seen_cores[record["match"]] = record["accession"]
            unique.append(record)
    cap = inherited["homolog_rows_max"]
    selected = core.farthest_first(receptor, unique, cap) if len(unique) > cap else unique
    selected_ids = {record["accession"] for record in selected}
    for record in unique:
        record["decision"] = "SELECTED" if record["accession"] in selected_ids else "EXCLUDE_DIVERSITY_CAP"
    write_csv("candidate_decisions.csv", records)
    write_csv("selected_homologs.csv", selected)
    save("full_alignment_evidence.json", alignment_evidence)
    receptor_a3m = ">Q92843_1_164\n" + receptor + "\n" + "".join(
        f">{r['accession']} taxon={r['taxon']} reviewed={r['reviewed']}\n{r['a3m']}\n" for r in selected)
    (OUT / "receptor_candidates.a3m").write_text(receptor_a3m, encoding="utf-8")
    (OUT / "selected_homologs_full.fasta").write_text("".join(
        f">{r['accession']} taxon={r['taxon']}\n{r['full_sequence']}\n" for r in selected), encoding="utf-8")
    content = core.complex_a3m(receptor, peptide, [(r["accession"], r["a3m"]) for r in selected])
    (OUT / "parent_candidate.a3m").write_text(content, encoding="utf-8")
    primary = core.diversity([receptor] + [r["match"] for r in selected])
    homolog_only = core.diversity([r["match"] for r in selected]) if selected else None
    gates = core.readiness(len(selected), primary, inherited["readiness_floor"])
    parser_checks = parse_with_installed_function(content, receptor, peptide, [r["a3m"] for r in selected])
    save("diversity_and_readiness.json", {"primary_query_inclusive": primary,
         "homolog_only_diagnostic": homolog_only, "frozen_floor": inherited["readiness_floor"],
         "gates": gates, "all_input_floor_gates_pass": all(gates.values()),
         "selected_row_order": ["Q92843_query"] + [r["accession"] for r in selected]})
    failed = [name for name, passed in gates.items() if not passed]
    ready = not failed and all(parser_checks.values())
    save("decision.json", {
        "audit_id": config["audit_id"], "status": "CPU_INPUT_AUDIT_COMPLETED",
        "source_rows": len(raw), "aligned_rows": len(alignment_evidence),
        "candidate_decision_counts": dict(Counter(r["decision"] for r in records)),
        "selected_nonquery_homologs": len(selected), "selected_reviewed": sum(r["reviewed"] for r in selected),
        "primary_clusters_80": primary["greedy_clusters_80"], "primary_neff_80": primary["neff_80"],
        "primary_minimum_column_nongap_fraction": primary["minimum_column_nongap_fraction"],
        "input_floor_gates": gates, "failed_input_gates": failed,
        "installed_parser_checks_pass": all(parser_checks.values()),
        "input_readiness": "GO_PREPARE_ONLY" if ready else "HOLD_INPUT_READINESS",
        "prediction_execution": "NOT_STARTED_NOT_AUTHORIZED_BY_THIS_AUDIT",
        "future_protocol_preparation_allowed": ready,
        "new_prediction_jobs": 0, "new_relaxation_jobs": 0, "new_reward_calls": 0,
        "new_downloads": 0, "new_training": 0,
        "interpretation": "Engineering input readiness only; no experimental affinity labels, no H/d validation and no demonstrated cause of prior low-confidence predictions.",
    })
    mismatches = [item["path"] for item in old_snapshot
                  if not Path(item["path"]).exists() or digest(Path(item["path"])) != item["sha256"]]
    save("preservation_after.json", {"files_checked": len(old_snapshot),
         "mismatches": mismatches, "all_match": not mismatches, "old_manifest_checks": old_manifests})
    assert not mismatches
    forbidden = sorted(name for name in sys.modules if name.split(".")[0] in {"torch", "transformers", "jax", "jaxlib", "colabfold", "esm", "openmm"})
    assert not forbidden
    output_bytes = sum(path.stat().st_size for path in OUT.rglob("*") if path.is_file())
    elapsed = time.monotonic() - started
    save("runtime_and_budget.json", {"started_utc": started_utc,
         "finished_utc": datetime.now(timezone.utc).isoformat(), "worker_wall_seconds": elapsed,
         "cpu_worker_wall_seconds_max": config["cpu_worker_wall_seconds_max"],
         "worker_wall_budget_pass": elapsed <= config["cpu_worker_wall_seconds_max"],
         "output_bytes_at_worker_completion": output_bytes, "output_bytes_max": config["output_bytes_max"],
         "output_budget_pass": output_bytes <= config["output_bytes_max"],
         "python": sys.version, "executable": sys.executable, "platform": platform.platform(),
         "biopython": Bio.__version__, "forbidden_runtime_modules_loaded": forbidden,
         "source_scope": "Existing local files only; no subprocess, network call, model inference, feature generation or GPU runtime import."})
    signal.alarm(0)
    print(json.dumps(load(OUT / "decision.json"), ensure_ascii=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        failure_path = OUT / "audit_failure.json"
        if not failure_path.exists():
            save("audit_failure.json", {"error_type": type(error).__name__, "message": str(error),
                 "new_predictions": 0, "new_reward_calls": 0, "status": "FAILURE_PRESERVED"})
        raise
