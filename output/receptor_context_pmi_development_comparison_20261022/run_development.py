"""Finite PMI development comparison; frozen inputs, one bounded GPU attempt."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OLD = ROOT / "output/receptor_context_mutation_probe_v0_20261022"
sys.path.insert(0, str(OLD))
from probe_core import (contrast, directory_bytes, mask_tokens, read_json,
                        run_bounded, sha256, write_json)
from probe_run import check_manifest, preservation_checks
from paired_core import (cache_key, paired_summary, parse_kd, parse_tables,
                         point_direction, predict_d, predict_h, sensitivity)

AA = "ACDEFGHIKLMNPQRSTVWY"
PRECISION_KEY = "FP32_eager_TF32off"
TOL = 1e-4
SOURCE = ROOT / "output/independent_affinity_validation_20261007/sources/PMC2856455.html"
SOURCE_TEXT = SOURCE.with_suffix(".txt")
CRITERIA = ROOT / "output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json"


def csv_write(name: str, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("Cannot infer CSV schema from zero rows")
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (HERE / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def preserve() -> list[dict]:
    return preservation_checks() + [check_manifest(OLD / "artifact_manifest.json")]


def code_paths() -> list[Path]:
    """Freeze the entrypoints and all imported local execution helpers."""
    return [HERE / name for name in ("plan_zh.md", "test_paired_core.py", "paired_core.py",
            "run_development.py", "deliver_results.py")] + [OLD / "probe_core.py", OLD / "probe_run.py"]


def audit_source() -> tuple[list[dict], list[dict]]:
    """Read every data row in local original tables; never inherit rounded ratios."""
    tables = parse_tables(SOURCE.read_text(encoding="utf-8-sig"))
    if len(tables) != 3:
        raise ValueError(f"Expected three original tables, found {len(tables)}")
    audit = []
    pairs = []
    for table_index, family, parent in ((0, "PMI", "TSFAEYWNLLSP"),
                                        (2, "P53", "ETFSDLWKLLPE")):
        source_rows = [row for row in tables[table_index] if len(row) >= 13]
        if not source_rows or "".join(source_rows[0][:12]).replace(" ", "") != parent:
            raise ValueError(f"Source WT sequence differs for {family}")
        parent_measure = parse_kd(source_rows[0][12])
        if parent_measure is None:
            raise ValueError("WT mean/SD missing")
        parent_kd, parent_sd = parent_measure
        for source_index, row in enumerate(source_rows, 1):
            sequence = "".join(row[:12]).replace(" ", "").replace("\xa0", "")
            measured = parse_kd(row[12]) if row[12].strip() else None
            item = {"source_table": table_index + 1, "source_data_row": source_index,
                    "family": family, "parent_sequence": parent, "variant_sequence": sequence,
                    "source_kd_sd_M_text": row[12], "kd_nM": measured[0] if measured else None,
                    "sd_nM": measured[1] if measured else None, "parent_kd_nM": parent_kd,
                    "parent_sd_nM": parent_sd, "ratio_from_means": measured[0] / parent_kd if measured else None,
                    "source_printed_ratio": row[13] if len(row) > 13 else "",
                    "point_direction": "not_defined", "site": None, "mutation_id": "",
                    "assay": "SPR competition", "receptor_construct": "MDM2(25-109)",
                    "termini": "free N amine / free C acid; explicitly neither N-acetylated nor C-amidated",
                    "formal_v2_label": "HOLD_NOT_RELEASED", "development_role": "EXPOSED_DEVELOPMENT",
                    "decision": ""}
            if not any(cell.strip() for cell in row):
                item["decision"] = "NONDATA_SEPARATOR"
            elif len(sequence) != 12 or any(aa not in AA for aa in sequence):
                item["decision"] = "EXCLUDE_LENGTH_OR_NONSTANDARD"
            elif sequence == parent:
                item["decision"] = "REFERENCE_WT_OR_SELF"
            else:
                changed = [i for i, (a, b) in enumerate(zip(parent, sequence), 1) if a != b]
                if len(changed) != 1:
                    item["decision"] = "EXCLUDE_NOT_SINGLE_SUBSTITUTION"
                else:
                    site = changed[0]
                    mutation = f"{parent[site-1]}{site if family == 'PMI' else site+16}{sequence[site-1]}"
                    item.update(site=site, mutation_id=f"{family}_{mutation}")
                    if family == "P53" and mutation == "P27A":
                        item["decision"] = "EXCLUDE_P27A_NO_RERUN"
                    elif measured is None:
                        item["decision"] = "HOLD_ND_NO_DOCUMENTED_DETECTION_LIMIT"
                    else:
                        direction = point_direction(measured[0], parent_kd)
                        item["point_direction"] = direction
                        if direction == "ambiguous":
                            item["decision"] = "EXCLUDE_BELOW_TWO_FOLD_POINT_EFFECT"
                        elif family == "P53":
                            item["decision"] = "HOLD_REPLICATE_COUNT_NOT_ATTRIBUTABLE_TO_P53"
                        else:
                            item["decision"] = "CONDITIONAL_DEVELOPMENT_ONLY"
                            item.update(sensitivity(*measured, parent_kd, parent_sd))
                            item["source_replicates"] = "at least three independent measurements; exact per-peptide n unavailable"
                            pairs.append(item.copy())
            audit.append(item)
    expected = {"PMI_S2A", "PMI_F3A", "PMI_E5A", "PMI_Y6A", "PMI_W7A", "PMI_N8A", "PMI_L10A"}
    if {pair["mutation_id"] for pair in pairs} != expected or len(pairs) != 7:
        raise ValueError(f"Unexpected source-derived PMI panel: {[p['mutation_id'] for p in pairs]}")
    if sum(pair["point_direction"] == "gain" for pair in pairs) != 1:
        raise ValueError("Expected one PMI point-effect gain")
    return audit, pairs


def verify_assets(config: dict) -> None:
    for asset in config["model_assets"]["assets"]:
        path = Path(asset["path"])
        if path.stat().st_size != asset["bytes"] or sha256(path) != asset["sha256"]:
            raise ValueError(f"Model asset hash/size mismatch: {path}")


def prepare() -> None:
    if (HERE / "run_config.json").exists():
        raise FileExistsError("Frozen config already exists; no implicit restart")
    tests = subprocess.run([sys.executable, "-B", str(HERE / "test_paired_core.py")], capture_output=True, text=True)
    (HERE / "test_before_run.log").write_text(tests.stdout + tests.stderr, encoding="utf-8")
    if tests.returncode:
        raise ValueError("Model-free scientific regressions failed")
    before = preserve()
    if any(item["mismatches"] for item in before):
        raise ValueError("Historical source/artifact integrity failed")
    if sha256(CRITERIA) != "23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b":
        raise ValueError("Formal v2 criteria changed")
    old_config = read_json(OLD / "run_config.json")
    if not read_json(OLD / "numeric_consistency.json")["passed"]:
        raise ValueError("Reused distributions did not pass prior numeric checks")
    if (old_config["precision"], old_config["attn_implementation"], old_config["tf32"],
            old_config["base_batch_size"], old_config["numeric_tolerance_nats"]) != ("FP32", "eager", False, 2, TOL):
        raise ValueError("Reuse precision/config differs")
    audit, pairs = audit_source()
    csv_write("source_input_audit.csv", audit)
    write_json(HERE / "pair_inventory.json", {"pairs": pairs, "primary_formal_metrics": None,
               "formal_label_release": False, "exposure": "DEVELOPMENT_CONDITIONAL_SENSITIVITY"})
    old_parents = {p["parent_id"]: p["full_sequence"] for p in old_config["inputs"]}
    old_contexts = {c["id"]: c["sequence"] for c in old_config["contexts"]}
    weight = next(a["sha256"] for a in old_config["assets"]["assets"] if a["path"].endswith("pytorch_model.bin"))
    cache = {}
    for line in (OLD / "raw_masked_distributions.jsonl").read_text().splitlines():
        record = json.loads(line)
        if record["phase"] != "base":
            continue
        key = cache_key(weight, old_contexts[record["context_id"]], old_parents[record["parent_id"]],
                        record["site"], PRECISION_KEY)
        if key in cache:
            raise ValueError("Duplicate context-specific cached base distribution")
        cache[key] = record
    jobs = []
    reused = []
    for pair in pairs:
        for context in old_config["contexts"]:
            job = {"job_id": f"base_{len(jobs)+1:03d}", "phase": "base", "site": pair["site"],
                   "parent_id": "PMI_WT_12MER", "context_id": context["id"], "mutation_id": pair["mutation_id"]}
            key = cache_key(weight, context["sequence"], pair["parent_sequence"], pair["site"], PRECISION_KEY)
            job["reuse"] = key in cache
            jobs.append(job)
            if key in cache:
                record = cache[key]
                if record["wt"] != pair["parent_sequence"][pair["site"]-1]:
                    raise ValueError("Cached WT token differs")
                reused.append({**record, **job, "reused_from_job_id": record["job_id"],
                               "full_cache_key": list(key)})
    targets = [job for job in jobs if job["context_id"] == old_config["contexts"][0]["id"]]
    for i, job in enumerate(reversed(targets), 1):
        jobs.append({**job, "job_id": f"repeat_{i:03d}", "phase": "repeat", "reuse": False})
    new_jobs = [job for job in jobs if not job["reuse"]]
    if (len(jobs), len(reused), len(new_jobs)) != (28, 9, 19):
        raise ValueError("Frozen 21-base / 9-reused / 19-new contract differs; no expansion")
    with (HERE / "reused_distributions.jsonl").open("x", encoding="utf-8") as stream:
        for record in reused:
            stream.write(json.dumps(record, allow_nan=False) + "\n")
    config = {
        "run_id": "PMI_CONTEXT_PAIRED_DEVELOPMENT_20261022",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User OK，开始吧 explicitly starts the prepared 7-pair/19-new-input development comparison.",
        "scope": "DEVELOPMENT_CONDITIONAL_SENSITIVITY; not formal v2 benchmark or independent validation",
        "parent_sequence": "TSFAEYWNLLSP", "parent_id": "PMI_WT_12MER", "contexts": old_config["contexts"],
        "pairs": pairs, "jobs": jobs, "new_jobs": new_jobs, "new_logical_inputs_max": 19,
        "reused_base_inputs": 9, "base_distributions_total": 21, "target_repeats": 7,
        "model_assets": old_config["assets"], "weight_sha256": weight,
        "precision_key": PRECISION_KEY, "fp32": True, "eval": True, "tf32": False, "autocast": False,
        "attn_implementation": "eager", "base_batch_size": 2, "repeat_batch_size": 1, "torch_seed": 20261022,
        "numeric_tolerance_nats": TOL, "formula": "d=logp(alt)-logp(WT); H=dt-(d1+d2)/2",
        "d_direction": "d>1e-4 gain; d<-1e-4 loss; otherwise tie",
        "h_direction": "both H_control>1e-4 gain; both <-1e-4 loss; otherwise abstain; all 7 stay in denominator",
        "raw_H_ablation": "signed H with +/-1e-4; report only, no result-dependent method selection",
        "old_gain_only_filter": "diagnostic only; never use to select the development evaluation denominator",
        "label_rule": "point ratio<=.5 gain or >=2 loss; sensitivity ranges conditional on n=3/t2/delta-method assumptions, not source CI",
        "formal_labels_released": 0, "primary_formal_metrics": None, "formal_gain_CI": None,
        "uncertainty_limit": "one exposed background and one point-effect gain; exact n/raw measurements/covariance absent",
        "gpu_process_wall_seconds_max": 120., "disk_bytes_max": 67108864,
        "new_structure_jobs": 0, "GraphPep_calls": 0, "FoldX_calls": 0, "ESM2_calls": 0, "training_steps": 0,
        "no_retry_or_precision_fallback": True,
        "preservation_before": before,
        "source_files": [{"path": str(path), "sha256": sha256(path)} for path in
                         (SOURCE, SOURCE_TEXT, CRITERIA, OLD / "run_config.json", OLD / "raw_masked_distributions.jsonl",
                          OLD / "artifact_manifest.json", HERE / "reused_distributions.jsonl", HERE / "pair_inventory.json")],
        "code_files": [{"path": str(path), "sha256": sha256(path)} for path in code_paths()],
        "gpu_before": subprocess.check_output(["nvidia-smi", "--query-gpu=name,memory.total,memory.used", "--format=csv,noheader"], text=True).strip(),
        "compute_apps_before": subprocess.check_output(["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader"], text=True).strip(),
    }
    verify_assets(config)
    write_json(HERE / "run_config.json", config)
    (HERE / "run_config.sha256").write_text(sha256(HERE / "run_config.json") + "  run_config.json\n")
    print(json.dumps({"status": "FROZEN_READY", "pairs": len(pairs), "reused": len(reused), "new_inputs": len(new_jobs),
                      "config_sha256": sha256(HERE / "run_config.json"), "gpu": config["gpu_before"]}), flush=True)


def validate() -> dict:
    config = read_json(HERE / "run_config.json")
    if sha256(HERE / "run_config.json") != (HERE / "run_config.sha256").read_text().split()[0]:
        raise ValueError("Frozen config changed")
    for row in config["source_files"] + config["code_files"]:
        if sha256(Path(row["path"])) != row["sha256"]:
            raise ValueError(f"Frozen file changed: {row['path']}")
    verify_assets(config)
    return config


def worker() -> None:
    config = validate()
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
    model_dir = config["model_assets"]["model_dir"]
    tokenizer = AutoTokenizer.from_pretrained(model_dir, local_files_only=True)
    aa_ids = {aa: tokenizer.convert_tokens_to_ids(aa) for aa in AA}
    if len(set(aa_ids.values())) != 20 or tokenizer.unk_token_id in aa_ids.values():
        raise ValueError("Canonical amino acid tokens invalid")
    contexts = {c["id"]: c["sequence"] for c in config["contexts"]}
    peptide = config["parent_sequence"]
    prepared = []
    for job in config["new_jobs"]:
        receptor = contexts[job["context_id"]]
        original = tokenizer.encode(receptor + peptide, add_special_tokens=True)
        tokens, index = mask_tokens(original, tokenizer.encode(receptor, add_special_tokens=False),
            tokenizer.encode(peptide, add_special_tokens=False), job["site"], tokenizer.cls_token_id,
            tokenizer.eos_token_id, tokenizer.mask_token_id)
        if len(tokens) != len(receptor) + len(peptide) + 2 or tokens.count(tokenizer.mask_token_id) != 1:
            raise ValueError("Bad token indexing")
        prepared.append((job, tokens, index))
    write_json(HERE / "token_index_checks.json", {"passed": len(prepared) == 19, "aa_token_ids": aa_ids,
        "indices": [{**job, "mask_token_index_0based": index} for job, _, index in prepared]})
    base = {(r["site"], r["context_id"]): r for r in
            [json.loads(line) for line in (HERE / "reused_distributions.jsonl").read_text().splitlines()]}
    print("Input checks passed; loading existing local FP32 PepMLM", flush=True)
    start = time.monotonic()
    model = AutoModelForMaskedLM.from_pretrained(model_dir, local_files_only=True,
                torch_dtype=torch.float32, attn_implementation="eager").eval().to("cuda")
    if any(p.is_floating_point() and p.dtype != torch.float32 for p in model.parameters()):
        raise ValueError("Model dtype differs from frozen FP32")
    runtime = {"torch": torch.__version__, "transformers": transformers.__version__, "cuda": torch.version.cuda,
               "parameter_dtype": "float32", "eval": not model.training, "autocast": False,
               "tf32_matmul": torch.backends.cuda.matmul.allow_tf32, "tf32_cudnn": torch.backends.cudnn.allow_tf32,
               "device": torch.cuda.get_device_name(), "model_load_wall_seconds": time.monotonic() - start,
               "logical_inputs_completed": 0}
    write_json(HERE / "runtime.json", runtime)
    count = 0
    checks = []
    with (HERE / "new_distributions.jsonl").open("x", encoding="utf-8") as stream:
        for phase, batch_size in (("base", 2), ("repeat", 1)):
            jobs = [item for item in prepared if item[0]["phase"] == phase]
            for offset in range(0, len(jobs), batch_size):
                batch = jobs[offset:offset+batch_size]
                if count + len(batch) > config["new_logical_inputs_max"]:
                    raise ValueError("New logical input cap exceeded")
                width = max(len(tokens) for _, tokens, _ in batch)
                ids = torch.full((len(batch), width), tokenizer.pad_token_id, dtype=torch.long, device="cuda")
                attention = torch.zeros_like(ids)
                for row, (_, tokens, _) in enumerate(batch):
                    ids[row, :len(tokens)] = torch.tensor(tokens, device="cuda")
                    attention[row, :len(tokens)] = 1
                write_json(HERE / "progress.json", {"attempted": count+len(batch), "completed": count,
                           "current_job_ids": [j["job_id"] for j, _, _ in batch]})
                with torch.inference_mode():
                    output = model(input_ids=ids, attention_mask=attention).logits
                    logits = torch.stack([output[row, index].float() for row, (_, _, index) in enumerate(batch)])
                    logp = logits.log_softmax(-1)
                if not torch.isfinite(logits).all() or not torch.isfinite(logp).all():
                    write_json(HERE / "nonfinite_distributions.json", {"job_ids": [j["job_id"] for j, _, _ in batch],
                        "logits_as_strings": [[str(v) for v in row] for row in logits.cpu().tolist()],
                        "logp_as_strings": [[str(v) for v in row] for row in logp.cpu().tolist()]})
                    raise ValueError("Nonfinite model output; no retry")
                logits_cpu, logp_cpu = logits.cpu().tolist(), logp.cpu().tolist()
                for row, (job, _, index) in enumerate(batch):
                    canonical = {aa: logp_cpu[row][aa_ids[aa]] for aa in AA}
                    wt = peptide[job["site"]-1]
                    record = {**job, "mask_token_index_0based": index, "wt": wt, "batch_size_actual": len(batch),
                              "canonical_log_probs": canonical, "logits_full_vocab": logits_cpu[row],
                              "log_probs_full_vocab": logp_cpu[row]}
                    stream.write(json.dumps(record, allow_nan=False) + "\n")
                    stream.flush()
                    count += 1
                    key = (job["site"], job["context_id"])
                    if phase == "base":
                        if key in base:
                            raise ValueError("Duplicate base distribution")
                        base[key] = record
                    else:
                        reference = base[key]["canonical_log_probs"]
                        delta_p = max(abs(canonical[aa] - reference[aa]) for aa in AA)
                        delta_d = max(abs((canonical[aa]-canonical[wt]) - (reference[aa]-reference[wt])) for aa in AA)
                        checks.append({"site": job["site"], "max_logp_delta": delta_p, "max_d_delta": delta_d,
                                       "reference_was_reused": base[key].get("reuse", False)})
                        if max(delta_p, delta_d) > TOL:
                            write_json(HERE / "numeric_failure.json", checks)
                            raise ValueError("Frozen numerical tolerance failed")
                write_json(HERE / "progress.json", {"attempted": count, "completed": count, "current_job_ids": []})
                print(f"Completed {count}/19 new logical inputs ({phase})", flush=True)
                del ids, attention, output, logits, logp
    runtime.update(logical_inputs_completed=count, peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated(),
                   peak_cuda_reserved_bytes=torch.cuda.max_memory_reserved())
    write_json(HERE / "runtime.json", runtime)
    write_json(HERE / "numeric_consistency.json", {"passed": len(checks) == 7 and count == 19,
               "tolerance_nats": TOL, "checks": checks,
               "max_logp_delta": max(c["max_logp_delta"] for c in checks),
               "max_d_delta": max(c["max_d_delta"] for c in checks)})
    print("All 19 new inputs and seven numerical repeats passed", flush=True)


def run() -> None:
    config = validate()
    if (HERE / "worker.log").exists() or (HERE / "budget_result.json").exists():
        raise FileExistsError("One model attempt already exists; no implicit retry")
    result = run_bounded([sys.executable, "-B", str(Path(__file__).resolve()), "worker"], HERE,
                          config["gpu_process_wall_seconds_max"], config["disk_bytes_max"], HERE / "worker.log")
    write_json(HERE / "budget_result.json", result)
    print(json.dumps(result), flush=True)
    if result["status"] != "completed":
        raise RuntimeError("Worker stopped; preserve partial evidence and HOLD")


def analyze() -> None:
    config = validate()
    budget = read_json(HERE / "budget_result.json")
    numeric = read_json(HERE / "numeric_consistency.json")
    if budget["status"] != "completed" or not numeric["passed"]:
        raise ValueError("No direction comparison after failed technical run")
    new = [json.loads(line) for line in (HERE / "new_distributions.jsonl").read_text().splitlines()]
    reused = [json.loads(line) for line in (HERE / "reused_distributions.jsonl").read_text().splitlines()]
    if len(new) != 19 or [r["job_id"] for r in new] != [j["job_id"] for j in config["new_jobs"]]:
        raise ValueError("Actual new inputs differ from frozen inventory")
    base_records = [r for r in new + reused if r["phase"] == "base"]
    base = {(r["site"], r["context_id"]): r["canonical_log_probs"] for r in base_records}
    if len(base) != 21 or len(base_records) != 21:
        raise ValueError("Base denominator differs")
    context_ids = [c["id"] for c in config["contexts"]]
    results = []
    max_self = 0.
    for pair in config["pairs"]:
        site = pair["site"]
        wt, alt = pair["parent_sequence"][site-1], pair["variant_sequence"][site-1]
        probs = [base[(site, cid)] for cid in context_ids]
        self_value = contrast(*probs, wt, wt)
        max_self = max(max_self, abs(self_value["h"]), abs(self_value["d_target"]))
        values = contrast(*probs, wt, alt)
        dp = predict_d(values["d_target"])
        hp = predict_h(values["h_control_1"], values["h_control_2"])
        raw_hp = predict_d(values["h"])
        results.append({"mutation_id": pair["mutation_id"], "site": site, "wt": wt, "alt": alt,
            "point_direction": pair["point_direction"], "ratio_from_means": pair["ratio"],
            "independent_sensitivity_lower": pair["independent_lower"], "independent_sensitivity_upper": pair["independent_upper"],
            "worst_correlation_sensitivity_lower": pair["worst_correlation_lower"],
            "worst_correlation_sensitivity_upper": pair["worst_correlation_upper"],
            "formal_label": "HOLD_NOT_RELEASED", **values, "d_prediction": dp, "h_prediction": hp,
            "raw_H_prediction_diagnostic_only": raw_hp, "d_matches_point_direction": dp == pair["point_direction"],
            "h_matches_point_direction": hp == pair["point_direction"],
            "raw_H_matches_point_direction": raw_hp == pair["point_direction"]})
    csv_write("paired_comparison.csv", results)
    comparison = paired_summary(results)
    raw_correct = sum(r["raw_H_matches_point_direction"] for r in results)
    after = preserve()
    if any(p["mismatches"] for p in after) or after != config["preservation_before"] or max_self > TOL:
        raise ValueError("Source preservation or self-substitution gate failed")
    summary = {
        "status": "DEVELOPMENT_CONDITIONAL_COMPARISON_COMPLETED",
        "development_pairs": 7, "point_effect_gains": 1, "point_effect_losses": 6,
        "new_logical_inputs": len(new), "reused_distributions": len(reused), "total_base_distributions": len(base),
        "numeric": numeric, "max_wt_self_difference": max_self, "comparison": comparison,
        "raw_H_sign_diagnostic_correct": raw_correct,
        "development_paired_gain_observed": comparison["h_correct"] > comparison["d_correct"],
        "decision": "EXPLORATORY_DEVELOPMENT_GAIN_ONLY" if comparison["h_correct"] > comparison["d_correct"] else "NO_OBSERVED_DEVELOPMENT_GAIN_HOLD_EXPANSION",
        "formal_labels_released": 0, "primary_formal_metrics": None, "formal_method_difference_95_CI": None,
        "paired_uncertainty": {"description": "Observed paired correctness table only; no credible formal gain CI from one exposed parent and one gain",
            "exact_n_and_raw_repeats": "not recovered", "independent_error_sensitivity_is_not_source_CI": True},
        "affinity_validation": "HOLD", "closed_loop": "HOLD", "training": 0, "new_structures": 0,
        "new_GraphPep_FoldX_ESM2_calls": 0, "budget": budget, "preservation_after": after,
        "claim_limit": "Direction agreement with seven exposed point-effect examples under stated uncertainty assumptions; not calibrated affinity or independent superiority"}
    write_json(HERE / "summary.json", summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("numeric", "preservation_after")},
                     ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("prepare", "run", "worker", "analyze"))
    mode = parser.parse_args().mode
    try:
        globals()[mode]()
    except Exception:
        if mode == "worker":
            write_json(HERE / "worker_failure.json", {"traceback": traceback.format_exc(),
                "partial_records_preserved": True, "no_retry": True})
        raise
