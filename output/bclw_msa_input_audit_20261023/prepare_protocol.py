"""Prepare exact future commands and inspect installed file I/O, without launch."""
import ast
import copy
import hashlib
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
INSTALLED = Path("/home/g/miniconda3/envs/colabfold_exp0/lib/python3.11/site-packages/colabfold")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path: Path):
    return json.loads(path.read_text())


def save(name: str, value) -> None:
    with (OUT / name).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2)
        stream.write("\n")


def main() -> None:
    decision = load(OUT / "decision.json")
    if decision["input_readiness"] != "GO_PREPARE_ONLY":
        raise RuntimeError("Input floor did not pass; do not prepare execution commands")
    prepared = OUT / "prepared_inputs/BCLW_Bim_A2eT_parent.a3m"
    prepared.parent.mkdir(exist_ok=True)
    with prepared.open("xb") as stream:
        stream.write((OUT / "parent_candidate.a3m").read_bytes())
    original = load(ROOT / "output/bclw_parent_kernel_ablation_v1/prepared_protocol.json")
    root_future = ROOT / "output/bclw_parent_local_msa_diagnostic_v1"
    assert not root_future.exists(), "Prepare a unique future output directory"
    correspondence = load(OUT / "input_correspondence.json")
    sequence_r, sequence_p = correspondence["receptor_sequence"], correspondence["peptide_sequence"]

    # Exercise exact installed I/O definitions while avoiding package/model imports.
    names = {"parse_fasta", "get_queries", "pair_sequences", "pad_sequences", "pair_msa", "msa_to_str"}
    source = (INSTALLED / "input.py").read_text()
    nodes = [node for node in ast.parse(source).body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert len(nodes) == len(names)
    snippet = "from __future__ import annotations\n\n" + "\n\n".join(ast.get_source_segment(source, node) for node in nodes)
    (OUT / "installed_input_used.py.txt").write_text(snippet, encoding="utf-8")
    namespace = {"Path": Path}
    exec(compile(snippet, str(INSTALLED / "input.py") + "::isolated_io", "exec"), namespace)
    queries, is_complex = namespace["get_queries"](prepared, sort_queries_by="none")
    parser_namespace = {"mk_mock_template": lambda sequence: {"mock_template_for": sequence}}
    exec(compile((OUT / "installed_parser_used.py.txt").read_text(), "isolated_batch_parser", "exec"), parser_namespace)
    unpaired, paired, query_seqs, cardinality, _ = parser_namespace["unserialize_msa"](queries[0][2], queries[0][1])
    roundtrip = namespace["msa_to_str"](unpaired, paired, query_seqs, cardinality)
    io_checks = {
        "exactly_one_parent_job": len(queries) == 1 and queries[0][0] == "BCLW_Bim_A2eT_parent",
        "query_sequence_exact": queries[0][1] == sequence_r + sequence_p,
        "custom_a3m_loaded_not_none": queries[0][2] == [prepared.read_text()],
        "complex_recognized": is_complex,
        "chain_sequences_preserved": query_seqs == [sequence_r, sequence_p],
        "installed_serialize_roundtrip_exact": roundtrip.strip() == prepared.read_text().strip(),
    }
    save("installed_io_validation.json", {
        "input_source_path": str(INSTALLED / "input.py"), "input_source_sha256": sha(INSTALLED / "input.py"),
        "checks": io_checks, "all_pass": all(io_checks.values()),
        "scope": "Installed get_queries/parse_fasta and A3M roundtrip only; no complete feature generation, JAX import or prediction was executed.",
    })
    assert all(io_checks.values()), io_checks
    commands = copy.deepcopy(original["commands"])
    for command in commands:
        seed = command["seed"]
        command["argv"][1] = str(prepared)
        command["argv"][2] = str(root_future / f"seed_{seed}/structures")
        index = command["argv"].index("--msa-mode")
        command["argv"][index + 1] = "mmseqs2_uniref_env"
        command["relax_argv"][-2] = str(root_future / f"seed_{seed}/structures")
        command["relax_argv"][-1] = str(root_future / f"seed_{seed}/relaxed")
    settings = copy.deepcopy(original["constant_prediction_settings"])
    settings["use_fast_kernels"] = False
    settings["msa_policy"] = {
        "input": "Local source-audited A3M containing 22 receptor homologs + receptor query; peptide query only.",
        "cli_mode": "mmseqs2_uniref_env",
        "enum_interpretation": "This installed CLI has no 'custom' enum; non-single mode preserves provided A3M. Exact .a3m file input bypasses get_msa_and_templates; the enum does NOT authorize MMseqs/remote search.",
        "remote_service_allowed": False, "templates_allowed": False,
        "requested_max_seq": 256, "requested_max_extra_seq": 512,
        "effective_difference": "Old single_sequence clamped these to 1/1 and discarded homolog rows. New non-single mode removes this clamp. This necessary MSA feature-depth change is part of the intervention, not a claim that only the file path changes.",
        "query_pairing": "Only the concatenated query is paired; no inferred receptor-peptide homolog pairing.",
    }
    baseline_quality = load(ROOT / "output/bclw_parent_kernel_ablation_v1/quality_gate_summary.json")
    baseline_rows = [{key: row[key] for key in ("seed", "iptm", "ptm", "receptor_mean_plddt", "peptide_mean_plddt", "interface_pae_max_directed_A", "all_diagnostic_gates_pass")} for row in baseline_quality["jobs"]]
    protocol = {
        "protocol_id": "BCLW_PARENT_LOCAL_MSA_DIAGNOSTIC_V1",
        "status": "PREPARED_NOT_STARTED_REQUIRES_FUTURE_EXPLICIT_START",
        "prepared_date_local": "2026-10-23 Asia/Shanghai",
        "question": "Under the completed no-fast-kernel baseline settings and original fixed seeds, is this small local receptor MSA sufficient to recover parent complex engineering readiness?",
        "not_tested": ["H/d effectiveness", "experimental affinity direction", "mutant effects", "full binding free energy", "general MSA utility or all causes of prior prediction failure"],
        "input_a3m": str(prepared), "input_a3m_sha256": sha(prepared),
        "input_representation": correspondence,
        "input_audit_decision_sha256": sha(OUT / "decision.json"),
        "inherited_msa_rules_sha256": sha(ROOT / "output/bclw_construct_msa_audit_20261022/scope_and_selection_frozen.json"),
        "implementation_freeze_sha256": sha(OUT / "selection_implementation_frozen.json"),
        "planned_output_root": str(root_future), "planned_output_root_absent_when_prepared": True,
        "baseline": {
            "directory": str(ROOT / "output/bclw_parent_kernel_ablation_v1"),
            "status": "COMPLETED_0_OF_2",
            "comparison": "Existing same-seed no-fast-kernel parent jobs; reuse confidence/geometry, never rerun or rescore baseline.",
            "rows": baseline_rows,
        },
        "intervention": "Local receptor MSA plus necessary non-single MSA ingestion/effective-depth changes; all other prediction and relaxation settings held fixed relative to kernel branch.",
        "prediction_settings": settings,
        "constant_relax_settings": original["constant_relax_settings"],
        "commands": commands,
        "code_and_asset_checks": {
            "batch_source_path": str(INSTALLED / "batch.py"), "batch_source_sha256": sha(INSTALLED / "batch.py"),
            "input_source_path": str(INSTALLED / "input.py"), "input_source_sha256": sha(INSTALLED / "input.py"),
            "relax_source_path": str(ROOT / "pipeline/relax_feedback_panel.py"),
            "relax_source_sha256": sha(ROOT / "pipeline/relax_feedback_panel.py"),
            "model_parameter_file": settings["model_parameter_file"],
            "model_parameter_expected_sha256": settings["model_parameter_sha256"],
            "model_parameter_verification": "Expected hash reused from prior actual run; bulk parameter file must be hashed again at future launch. No weight loading/hash claimed in this input audit.",
        },
        "budget": {
            "new_prediction_jobs_max": 2, "raw_final_structures_max": 2,
            "restrained_relaxed_structures_max": 2, "run_serially": True,
            "gpu_tagged_process_wall_seconds_max": 1800,
            "cpu_relax_wall_seconds_max": 900, "disk_bytes_max": 536870912,
            "accounting": "Sum each child wall time once; wall time is not GPU-utilization integral. Watchdog enforces aggregate active+completed time and output cap.",
            "new_reward_calls_max": 0, "new_downloads_max": 0, "training_jobs_max": 0,
        },
        "quality_gate_unchanged": original["quality_gate_unchanged"],
        "other_required_gates": original["other_required_gates"],
        "additional_diagnostic_gate": original["additional_diagnostic_gate"],
        "analysis_freeze": {
            "chain_mapping": "A=unchanged 164-aa receptor; B=unchanged 23-aa peptide; exact sequences and contiguous 1-based numbering in raw and relaxed PDB.",
            "confidence": "Raw score JSON ipTM; raw CA pLDDT means for A/B; maximum over ALL A->B and B->A PAE entries, not contact-restricted or median.",
            "geometry": "After identical CPU restrained relaxation, minimum cross-chain heavy-atom distance >=1.5 A and >=1 unique cross-chain residue pair with any heavy atoms <=6 A. Also retain raw geometry and heavy-atom pair count separately.",
            "consistency": "Cross-check JSON pLDDT against raw CA B factors and score PAE against exported PAE; report exact tolerances from serialization, fail on inconsistent identity/array shapes/nonfinite values.",
            "reference_code": str(ROOT / "output/bclw_parent_kernel_ablation_v1/summarize_kernel_ablation.py"),
            "reference_code_sha256": sha(ROOT / "output/bclw_parent_kernel_ablation_v1/summarize_kernel_ablation.py"),
            "reference_use": "Reuse metric definitions or isolated pure helpers in a NEW output-directed evaluator. Existing runner and summarizer contain old hardcoded output paths and must not be directly launched for this protocol.",
            "paired_changes": "For each fixed seed report new-minus-existing-baseline confidence and geometry; never pool two seeds into a pass.",
        },
        "preflight_before_future_start": [
            "Require a future explicit instruction to start this prepared protocol. Current authorization is preparation only.",
            "Verify applicable AGENTS, protocol/input SHA256, old source/pilot/kernel integrity and inherited unchanged frozen criteria.",
            "Verify exact installed code/version/executable/model hashes, GPU availability and disk capacity; do not reinstall, download, kill unrelated jobs or change precision/options on OOM.",
            "Require a clean new output root. Freeze any necessary unique path change before launch and preserve old outputs.",
            "Repeat installed A3M ingestion checks and ensure a3m_lines is provided and templates disabled. Deny any fallback to remote MSA/template services or weight downloads.",
            "Use a new budget wrapper/evaluator bound to this protocol and new output root. Do not directly reuse old scripts with hardcoded kernel paths. Count each child interval once.",
            "Verify saved A3M still contains 22 receptor homologs, query-only peptide, unchanged chain sequences, and no single_sequence clamp. Record actual config and MSA feature depth; lost homologs are an execution failure.",
        ],
        "frozen_decision": {
            "2_of_2_all_original_and_auxiliary_gates_pass": "GO_TO_SEPARATE_REVIEW_ONLY; evidence of parent readiness under this MSA+settings bundle, not affinity or H validation. No mutant run automatically follows.",
            "1_of_2_pass": "HOLD_SEED_SENSITIVE; keep denominator 2, no reroll.",
            "0_of_2_pass": "HOLD_THIS_MSA_DIAGNOSTIC; no gate loosening or expansion; does not prove that all other MSAs, constructs or runtime factors are ineffective.",
            "input_runtime_budget_or_execution_failure": "Stop, preserve failed jobs/logs and denominator; no automatic fallback, new seed, precision, template or additional compute.",
            "better_confidence_but_gate_failure": "Descriptive improvement only, classified as gate failure.",
        },
        "scope_exclusions": original["excluded"] + ["remote MSA service", "template search", "baseline rerun"],
        "formal_direction_criteria": "Unchanged v2. No formal Kd labels or direction metrics are released. Existing failures/ties and GraphPep/FoldX 0.05 bands remain preserved.",
    }
    save("prepared_diagnostic_protocol.json", protocol)
    (OUT / "prepared_diagnostic_protocol.sha256").write_text(sha(OUT / "prepared_diagnostic_protocol.json") + "\n")
    save("prepared_protocol_checks.json", {
        "input_copied_exactly": sha(prepared) == sha(OUT / "parent_candidate.a3m"),
        "all_executable_paths_exist": all(Path(command["argv"][0]).is_file() and Path(command["relax_argv"][0]).is_file() for command in commands),
        "commands_prepared": len(commands), "commands_launched": 0,
        "new_output_directory_created": root_future.exists(),
        "protocol_sha256": sha(OUT / "prepared_diagnostic_protocol.json"),
    })
    print(json.dumps({"protocol": protocol["protocol_id"], "status": protocol["status"], "io_checks": io_checks, "commands_launched": 0}, ensure_ascii=False))


if __name__ == "__main__":
    main()
