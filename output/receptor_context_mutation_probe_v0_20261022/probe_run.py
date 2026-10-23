"""Bounded PepMLM context probe. No structure/reward calls or sampling."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import time
import traceback

from probe_core import (build_jobs, contrast, directory_bytes, mask_tokens, read_json,
                        run_bounded, select_suggestions, sha256, write_json)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PLAN = ROOT / "output/algorithm_hypothesis_plan_20261022"
AA = "ACDEFGHIKLMNPQRSTVWY"
TOL = 1e-4


def check_manifest(manifest: Path) -> dict:
    failures = []
    data = read_json(manifest)
    for row in data["files"]:
        path = Path(row["path"])
        if not path.is_absolute():
            path = manifest.parent / path
        if not path.exists() or sha256(path) != row["sha256"]:
            failures.append(str(path))
    return {"manifest": str(manifest), "manifest_sha256": sha256(manifest),
            "files": len(data["files"]), "mismatches": failures}


def preservation_checks() -> list[dict]:
    return [check_manifest(path) for path in (
        PLAN / "artifact_manifest.json",
        ROOT / "output/bh3_bclw_six_complex_pilot_20261020/run_result_manifest_v1.json",
        ROOT / "output/bclw_parent_kernel_ablation_v1/artifact_manifest.json")]


def write_csv(name: str, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("Explicit columns required for empty CSV")
    with (HERE / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def prepare() -> None:
    """Hash and freeze before loading weights; never overwrite a frozen run."""
    config_path = HERE / "run_config.json"
    if config_path.exists():
        raise FileExistsError("Run config is already frozen; no automatic restart")
    tests = subprocess.run([sys.executable, "-B", str(HERE / "test_probe.py")],
                           capture_output=True, text=True)
    (HERE / "test_before_run.log").write_text(tests.stdout + tests.stderr, encoding="utf-8")
    if tests.returncode:
        raise ValueError("Model-free regressions failed before freeze")
    preservation = preservation_checks()
    if any(item["mismatches"] for item in preservation):
        raise ValueError(f"Source manifests changed: {preservation}")
    assets = read_json(PLAN / "existing_model_asset_manifest.json")
    for item in assets["assets"]:
        path = Path(item["path"])
        if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise ValueError(f"Model hash mismatch: {path}")
    reused = read_json(PLAN / "source_reuse_manifest.json")
    for item in reused["files"]:
        if sha256(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"Reused source changed: {item['path']}")
    inputs = read_json(PLAN / "development_input_inventory.json")["inputs"]
    contexts = read_json(PLAN / "context_controls.json")["controls"]
    fasta = "".join(f">{item['id']}\n{item['sequence']}\n" for item in contexts)
    if (PLAN / "contexts_prepared.fasta").read_text() != fasta:
        raise ValueError("Prepared FASTA and JSON differ")
    if len({item["sequence"] for item in contexts}) != 3:
        raise ValueError("Controls must be distinct")
    for item in contexts:
        if Counter(item["sequence"]) != Counter(contexts[0]["sequence"]):
            raise ValueError("Control composition differs")
    for parent in inputs:
        if parent["receptor_sequence_for_PepMLM"] != contexts[0]["sequence"]:
            raise ValueError("Receptor input differs")
        if set(parent["eligible_positions"]) & set(parent["protected_positions"] +
                parent["unresolved_positions"] + parent["protocol_excluded_positions"]):
            raise ValueError("Forbidden mask site")
    jobs = build_jobs(inputs, [item["id"] for item in contexts])
    if len(jobs) != 60:
        raise ValueError("Frozen logical input count differs")
    baseline_paths = [ROOT / "output/feedback_panel_mdm2_20261005/native_graphpep/graphpep_inputs" /
                      sub / "residue_feedback.csv" for sub in
                      ("0000_P53WT_1YCR_1YCR", "0001_PMI_3EQS_3EQS")]
    baseline = [{"parent_id": parent["parent_id"], "path": str(path), "sha256": sha256(path)}
                for parent, path in zip(inputs, baseline_paths)]
    command = ["nvidia-smi", "--query-gpu=name,memory.total,memory.used", "--format=csv,noheader"]
    gpu = subprocess.check_output(command, text=True)
    active = subprocess.check_output(["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory",
                                      "--format=csv,noheader"], text=True)
    config = {
        "run_id": "RECEPTOR_CONTEXT_MUTATION_PROBE_V0_20261022",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User OK试试吧 authorizes this prepared <=60-input technical experiment only.",
        "model_free_test_returncode": tests.returncode,
        "proposal": str(PLAN / "proposal.json"), "proposal_sha256": sha256(PLAN / "proposal.json"),
        "historical_proposal_status": "NOT_RUN refers to preparation time; execution is recorded here.",
        "input_files": [{"path": str(PLAN / name), "sha256": sha256(PLAN / name)} for name in
                        ("development_input_inventory.json", "context_controls.json", "contexts_prepared.fasta")],
        "inputs": inputs, "contexts": contexts, "jobs": jobs, "assets": assets,
        "precision": "FP32", "tf32": False, "autocast": False, "training": False,
        "attn_implementation": "eager", "base_batch_size": 2, "repeat_batch_size": 1,
        "repeat_order": "all target jobs reversed", "numeric_tolerance_nats": TOL,
        "numeric_gate": "maximum canonical-amino-acid masked log-prob and derived d differences <=1e-4",
        "torch_seed": 20261022, "random_baseline_seed": 2026102203,
        "candidate_pool": "19 non-WT canonical substitutions per eligible site; remove named exposed mutants below",
        "historical_candidate_exclusions": [
            {"parent_id": "P53WT_1YCR", "site": 6, "wt": "L", "alt": "A", "name": "P53_L22A"},
            {"parent_id": "PMI_3EQS", "site": 8, "wt": "N", "alt": "A", "name": "PMI_N8A"}],
        "other_historical_exclusions": "P53_P27A site 11 and PMI_W7A site 7 are not masked at all.",
        "formula": "d=logp(alt)-logp(WT); H=d_target-(d_shuffle1+d_shuffle2)/2",
        "filter": "d_target>1e-4 and H_control_1>1e-4 and H_control_2>1e-4",
        "control_disagreement": "H_control_1 and H_control_2 have opposite signs beyond +/-1e-4",
        "ranking": "H descending, site ascending, alt alphabetical; <=1 per site, <=6 per parent",
        "baseline_ranking": {
            "unadjusted_d": "same 283 candidate pool; d_target descending, site, alt; no H-based filter",
            "static_GraphPep": "existing NATIVE-pose residue means ascending then site, alt alphabetical; position baseline only, no amino-acid preference",
            "uniform_random": "shuffle exact same pool per parent with Python Random(2026102203); no H-based filter",
            "top_k_comparison": "for each parent k equals H suggestion count; select <=1 per site for all arms; report exact-mutant and site overlap descriptively"},
        "baseline_sources": baseline, "cache": "none; each logical input executes independently",
        "gpu_process_wall_seconds_max": 600., "disk_bytes_max": 67108864,
        "logical_inputs_max": 60, "new_structures": 0, "reward_calls": 0,
        "claim_limit": "Technical context sensitivity only; no affinity, specificity, novelty or superiority validation.",
        "preservation_before": preservation,
        "source_reuse_manifest": reused,
        "gpu_before": gpu.strip(), "compute_apps_before": active.strip(),
        "code_hashes": [{"path": str(HERE / name), "sha256": sha256(HERE / name)}
                        for name in ("probe_core.py", "probe_run.py", "test_probe.py")],
    }
    write_json(config_path, config)
    (HERE / "run_config.sha256").write_text(sha256(config_path) + "  run_config.json\n")
    print(json.dumps({"status": "FROZEN_READY", "jobs": len(jobs), "config_sha256": sha256(config_path),
                      "preservation": preservation, "gpu": gpu.strip()}, ensure_ascii=False), flush=True)


def validate_config() -> dict:
    config = read_json(HERE / "run_config.json")
    if sha256(HERE / "run_config.json") != (HERE / "run_config.sha256").read_text().split()[0]:
        raise ValueError("Frozen run config changed")
    for item in config["input_files"] + config["code_hashes"] + config["baseline_sources"]:
        if sha256(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"Frozen file changed: {item['path']}")
    if sha256(Path(config["proposal"])) != config["proposal_sha256"]:
        raise ValueError("Frozen proposal changed")
    for item in config["assets"]["assets"]:
        path = Path(item["path"])
        if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
            raise ValueError(f"Frozen model asset changed: {path}")
    return config


def worker() -> None:
    """Only mode that loads the model; guarded externally by run_bounded."""
    config = validate_config()
    os.environ.update(HF_HUB_OFFLINE="1", TRANSFORMERS_OFFLINE="1", TOKENIZERS_PARALLELISM="false",
                      CUBLAS_WORKSPACE_CONFIG=":4096:8")
    import torch
    import transformers
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    torch.set_num_threads(2)
    torch.manual_seed(config["torch_seed"])
    torch.cuda.manual_seed_all(config["torch_seed"])
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    torch.use_deterministic_algorithms(True)
    model_dir = config["assets"]["model_dir"]
    tokenizer = AutoTokenizer.from_pretrained(model_dir, local_files_only=True)
    aa_ids = {aa: tokenizer.convert_tokens_to_ids(aa) for aa in AA}
    if len(set(aa_ids.values())) != 20 or tokenizer.unk_token_id in aa_ids.values():
        raise ValueError("Canonical token mapping invalid")
    parents = {p["parent_id"]: p for p in config["inputs"]}
    contexts = {c["id"]: c["sequence"] for c in config["contexts"]}
    prepared = []
    for job in config["jobs"]:
        peptide = parents[job["parent_id"]]["full_sequence"]
        receptor = contexts[job["context_id"]]
        ids = tokenizer.encode(receptor + peptide, add_special_tokens=True)
        rids = tokenizer.encode(receptor, add_special_tokens=False)
        pids = tokenizer.encode(peptide, add_special_tokens=False)
        masked, index = mask_tokens(ids, rids, pids, job["site"], tokenizer.cls_token_id,
                                    tokenizer.eos_token_id, tokenizer.mask_token_id)
        if len(ids) != len(receptor) + len(peptide) + 2 or masked.count(tokenizer.mask_token_id) != 1:
            raise ValueError("Tokenization is not one residue per token")
        prepared.append((job, masked, index))
    if len(prepared) != 60:
        raise ValueError("Input count differs before loading weights")
    write_json(HERE / "token_index_checks.json", {"checked": len(prepared), "passed": True,
               "aa_token_ids": aa_ids, "indices": [{**j, "mask_token_index_0based": i,
                "nonpad_token_count": len(ids)} for j, ids, i in prepared]})
    print("Token/index checks passed for 60 logical inputs; loading local FP32 model", flush=True)
    started = time.monotonic()
    model = AutoModelForMaskedLM.from_pretrained(model_dir, local_files_only=True,
                torch_dtype=torch.float32, attn_implementation="eager").eval().to("cuda")
    if any(p.is_floating_point() and p.dtype != torch.float32 for p in model.parameters()):
        raise ValueError("Model parameter precision differs from frozen FP32")
    runtime = {"torch": torch.__version__, "transformers": transformers.__version__,
               "cuda": torch.version.cuda, "device": torch.cuda.get_device_name(),
               "parameter_dtype": "float32", "eval": not model.training, "autocast": False,
               "tf32_matmul": torch.backends.cuda.matmul.allow_tf32,
               "tf32_cudnn": torch.backends.cudnn.allow_tf32, "attention": "eager",
               "model_load_wall_seconds": time.monotonic() - started, "logical_inputs_completed": 0}
    write_json(HERE / "runtime.json", runtime)
    target_base = {}
    repeat_diffs = []
    count = 0
    with (HERE / "raw_masked_distributions.jsonl").open("x", encoding="utf-8") as stream:
        for phase, batch_size in (("base", 2), ("repeat", 1)):
            phase_jobs = [item for item in prepared if item[0]["phase"] == phase]
            for offset in range(0, len(phase_jobs), batch_size):
                batch = phase_jobs[offset:offset + batch_size]
                if count + len(batch) > config["logical_inputs_max"]:
                    raise ValueError("Logical input cap exceeded")
                width = max(len(ids) for _, ids, _ in batch)
                ids_tensor = torch.full((len(batch), width), tokenizer.pad_token_id,
                                        dtype=torch.long, device="cuda")
                attention = torch.zeros_like(ids_tensor)
                for row, (_, ids, _) in enumerate(batch):
                    ids_tensor[row, :len(ids)] = torch.tensor(ids, dtype=torch.long, device="cuda")
                    attention[row, :len(ids)] = 1
                write_json(HERE / "progress.json", {"logical_inputs_attempted": count + len(batch),
                           "logical_inputs_completed": count, "current_job_ids": [j["job_id"] for j, _, _ in batch]})
                with torch.inference_mode():
                    output = model(input_ids=ids_tensor, attention_mask=attention).logits
                    logits = torch.stack([output[row, index, :].float()
                                          for row, (_, _, index) in enumerate(batch)])
                    log_probs = logits.log_softmax(dim=-1)
                if not torch.isfinite(logits).all() or not torch.isfinite(log_probs).all():
                    write_json(HERE / "nonfinite_distributions.json", {
                        "job_ids": [j["job_id"] for j, _, _ in batch],
                        "logits_as_strings": [[str(v) for v in row] for row in logits.cpu().tolist()],
                        "log_probs_as_strings": [[str(v) for v in row] for row in log_probs.cpu().tolist()],
                        "failed_logical_inputs": len(batch), "completed_finite_records": count})
                    raise ValueError("Nonfinite model distribution")
                raw_logits = logits.cpu().tolist()
                raw_logp = log_probs.cpu().tolist()
                for row, (job, _, index) in enumerate(batch):
                    peptide = parents[job["parent_id"]]["full_sequence"]
                    canonical = {aa: raw_logp[row][aa_ids[aa]] for aa in AA}
                    key = (job["parent_id"], job["site"])
                    record = {**job, "mask_token_index_0based": index, "wt": peptide[job["site"] - 1],
                              "batch_size_actual": len(batch), "logits_full_vocab": raw_logits[row],
                              "log_probs_full_vocab": raw_logp[row], "canonical_log_probs": canonical}
                    stream.write(json.dumps(record, allow_nan=False) + "\n")
                    stream.flush()
                    count += 1
                    if phase == "base" and job["context_id"] == config["contexts"][0]["id"]:
                        target_base[key] = record
                    if phase == "repeat":
                        prior = target_base[key]
                        delta_p = max(abs(canonical[a] - prior["canonical_log_probs"][a]) for a in AA)
                        wt = record["wt"]
                        delta_d = max(abs((canonical[a] - canonical[wt]) -
                                          (prior["canonical_log_probs"][a] - prior["canonical_log_probs"][wt]))
                                      for a in AA)
                        repeat_diffs.append({"parent_id": key[0], "site": key[1],
                                             "max_logp_delta": delta_p, "max_d_delta": delta_d})
                        if max(delta_p, delta_d) > TOL:
                            write_json(HERE / "numeric_failure.json", repeat_diffs)
                            raise ValueError("Frozen batch/order numerical tolerance failed; no rescue")
                print(f"Completed {count}/60 logical inputs ({phase})", flush=True)
                write_json(HERE / "progress.json", {"logical_inputs_attempted": count,
                           "logical_inputs_completed": count, "current_job_ids": []})
                del output, logits, log_probs, ids_tensor, attention
    runtime.update(logical_inputs_completed=count, peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated(),
                   peak_cuda_reserved_bytes=torch.cuda.max_memory_reserved(),
                   post_load_wall_seconds=time.monotonic() - started - runtime["model_load_wall_seconds"])
    write_json(HERE / "runtime.json", runtime)
    write_json(HERE / "numeric_consistency.json", {"passed": len(repeat_diffs) == 15,
               "tolerance_nats": TOL, "checks": repeat_diffs,
               "max_logp_delta": max(x["max_logp_delta"] for x in repeat_diffs),
               "max_d_delta": max(x["max_d_delta"] for x in repeat_diffs)})
    print("All 60 inputs completed; numerical consistency passed", flush=True)


def run() -> None:
    config = validate_config()
    if (HERE / "budget_result.json").exists() or (HERE / "worker.log").exists():
        raise FileExistsError("A model attempt already exists; no automatic rerun")
    result = run_bounded([sys.executable, "-B", str(Path(__file__).resolve()), "worker"], HERE,
                          config["gpu_process_wall_seconds_max"], config["disk_bytes_max"], HERE / "worker.log")
    write_json(HERE / "budget_result.json", result)
    print(json.dumps(result), flush=True)
    if result["status"] != "completed":
        raise RuntimeError("Worker stopped; preserve partial raw evidence, do not rerun")


def analyze() -> None:
    config = validate_config()
    budget = read_json(HERE / "budget_result.json")
    numeric = read_json(HERE / "numeric_consistency.json")
    if budget["status"] != "completed" or not numeric["passed"]:
        raise ValueError("No suggestions allowed after a failed run")
    records = [json.loads(line) for line in (HERE / "raw_masked_distributions.jsonl").read_text().splitlines()]
    if len(records) != 60 or [r["job_id"] for r in records] != [j["job_id"] for j in config["jobs"]]:
        raise ValueError("Raw job count/order differs")
    base = {(r["parent_id"], r["site"], r["context_id"]): r for r in records if r["phase"] == "base"}
    excluded = {(r["parent_id"], r["site"], r["alt"]) for r in config["historical_candidate_exclusions"]}
    context_ids = [c["id"] for c in config["contexts"]]
    rows = []
    max_self = 0.
    for parent in config["inputs"]:
        for site in parent["eligible_positions"]:
            wt = parent["full_sequence"][site - 1]
            probs = [base[(parent["parent_id"], site, cid)]["canonical_log_probs"] for cid in context_ids]
            zero = contrast(*probs, wt, wt)
            max_self = max(max_self, abs(zero["h"]), abs(zero["d_target"]))
            for alt in AA:
                if alt == wt:
                    continue
                values = contrast(*probs, wt, alt)
                rows.append({"parent_id": parent["parent_id"], "site": site, "wt": wt, "alt": alt,
                             **values, "excluded_historical": (parent["parent_id"], site, alt) in excluded})
    write_csv("all_substitution_diagnostics.csv", rows)
    pool = [r for r in rows if not r["excluded_historical"]]
    suggestions = select_suggestions(pool)
    fields = ["parent_id", "site", "wt", "alt", "d_target", "d_control_1", "d_control_2",
              "h_control_1", "h_control_2", "h", "passes_filter", "control_sign_disagreement", "excluded_historical"]
    with (HERE / "technical_suggestions.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: r[k] for k in fields} for r in suggestions)
    rng = random.Random(config["random_baseline_seed"])
    baseline_rows = []
    parents_summary = []
    for parent in config["inputs"]:
        pid = parent["parent_id"]
        group = [r for r in pool if r["parent_id"] == pid]
        selected = [r for r in suggestions if r["parent_id"] == pid]
        source = next(s for s in config["baseline_sources"] if s["parent_id"] == pid)
        with Path(source["path"]).open(newline="") as stream:
            profile = {int(r["peptide_position"]): float(r["mean"])
                       for r in csv.DictReader(stream) if int(r["count"]) > 0}
        orders = {
            "unadjusted_d": sorted(group, key=lambda r: (-r["d_target"], r["site"], r["alt"])),
            "static_GraphPep_native_position_only": sorted(group, key=lambda r: (profile[r["site"]], r["site"], r["alt"])),
            "uniform_random": group.copy()}
        rng.shuffle(orders["uniform_random"])
        selected_keys = {(r["site"], r["alt"]) for r in selected}
        overlaps = {}
        for arm, ordered in orders.items():
            seen = set()
            chosen = []
            for row in ordered:
                if row["site"] not in seen and len(chosen) < len(selected):
                    chosen.append(row)
                    seen.add(row["site"])
            exact = len(selected_keys & {(r["site"], r["alt"]) for r in chosen})
            site_overlap = len({r["site"] for r in selected} & {r["site"] for r in chosen})
            overlaps[arm] = {"k": len(chosen), "exact_mutant_overlap":
                             None if arm == "static_GraphPep_native_position_only" else exact,
                             "site_overlap": site_overlap}
            position_only = arm == "static_GraphPep_native_position_only"
            baseline_rows.extend({"parent_id": pid, "arm": arm, "rank": rank,
                                  "site": r["site"], "wt": r["wt"],
                                  "alt": "" if position_only else r["alt"],
                                  "d_target": "" if position_only else r["d_target"],
                                  "h": "" if position_only else r["h"]}
                                 for rank, r in enumerate(chosen, 1))
        parents_summary.append({"parent_id": pid, "mask_sites": len(parent["eligible_positions"]),
            "candidate_pool": len(group), "filter_pass_candidates": sum(r["passes_filter"] for r in group),
            "filter_pass_sites": len({r["site"] for r in group if r["passes_filter"]}),
            "target_positive_candidates": sum(r["d_target"] > TOL for r in group),
            "control_sign_disagreement": sum(r["control_sign_disagreement"] for r in group),
            "suggestions": len(selected), "overlaps": overlaps,
            "h_min": min(r["h"] for r in group), "h_max": max(r["h"] for r in group),
            "h_control_1_abs_max": max(abs(r["h_control_1"]) for r in group),
            "h_control_2_abs_max": max(abs(r["h_control_2"]) for r in group)})
    if baseline_rows:
        write_csv("baseline_top_k_diagnostics.csv", baseline_rows)
    preserved = preservation_checks()
    if any(item["mismatches"] for item in preserved) or max_self > TOL:
        raise ValueError("Integrity/self substitution verification failed")
    context_effect = max(abs(r[k]) for r in pool for k in ("h_control_1", "h_control_2"))
    summary = {"status": "TECHNICAL_PROBE_COMPLETED", "logical_inputs": len(records),
        "base_inputs": 45, "repeat_inputs": 15, "substitutions_before_exclusions": len(rows),
        "historical_exclusions": len(rows) - len(pool), "candidate_pool": len(pool),
        "max_wt_self_difference": max_self, "numeric": numeric,
        "context_effect_abs_max_nats": context_effect, "context_order_effect_above_numeric_tolerance": context_effect > TOL,
        "parents": parents_summary, "budget": budget,
        "decision": "GO_TECHNICAL_DEVELOPMENT_ONLY" if context_effect > TOL and suggestions else "HOLD_EXPANSION",
        "affinity_validation": "HOLD: qualified independent bidirectional labels/protocols not met",
        "closed_loop": "HOLD: no child sequence generation or structural evaluation authorized here",
        "preservation_after": preserved, "claim_limit": config["claim_limit"]}
    write_json(HERE / "summary.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("numeric", "preservation_after")},
                     ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "run", "worker", "analyze"))
    args = parser.parse_args()
    try:
        globals()[args.mode]()
    except Exception:
        if args.mode == "worker":
            write_json(HERE / "worker_failure.json", {"traceback": traceback.format_exc(),
                "raw_records_preserved": True, "no_automatic_retry": True})
        raise
