"""CPU-only audit of existing distributions and bounded local experimental pairs."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PROBE = ROOT / "output/receptor_context_mutation_probe_v0_20261022"
PAIRED = ROOT / "output/receptor_context_pmi_development_comparison_20261022"
SOURCES = ROOT / "output/independent_affinity_validation_20261007"
SCREEN = ROOT / "output/reward_direction_source_screen_20261020"
sys.path.insert(0, str(PROBE))
from probe_core import contrast, directory_bytes, read_json, run_bounded, sha256, write_json
from probe_run import check_manifest, preservation_checks
sys.path.insert(0, str(PAIRED))
from paired_core import predict_d, predict_h
from audit_core import disagreement_type, distribution_key, exact_pair, near_pair_reason, variant_sequence

AA = "ACDEFGHIKLMNPQRSTVWY"
PMI = "TSFAEYWNLLSP"
PMI_L10A = "TSFAEYWNLASP"
DISAGREEMENT_TYPES = {"strict_opposite", "H_abstains", "d_tie_H_directional", "method_failure"}


def read_csv(path: Path, delimiter: str = ",") -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(name: str, rows: list[dict], empty_fields: list[str] | None = None) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row)) if rows else empty_fields
    if not fields:
        raise ValueError("An empty table must retain an explicit schema")
    with (HERE / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def history() -> list[dict]:
    return preservation_checks() + [check_manifest(PROBE / "artifact_manifest.json"),
                                   check_manifest(PAIRED / "artifact_manifest.json")]


def code_paths() -> list[Path]:
    return [HERE / name for name in ("plan_zh.md", "test_audit.py", "audit_core.py", "audit_disagreements.py")] + [
        PROBE / "probe_core.py", PROBE / "probe_run.py", PAIRED / "paired_core.py"]


def freeze() -> None:
    if (HERE / "audit_config.json").exists():
        raise FileExistsError("Audit is already frozen")
    tests = subprocess.run([sys.executable, "-B", str(HERE / "test_audit.py")], capture_output=True, text=True)
    (HERE / "test_before_run.log").write_text(tests.stdout + tests.stderr, encoding="utf-8")
    if tests.returncode:
        raise ValueError("Identity/direction tests failed")
    preserved = history()
    if any(item["mismatches"] for item in preserved):
        raise ValueError("Prior artifacts do not match their original manifests")
    input_paths = [PROBE / name for name in ("run_config.json", "raw_masked_distributions.jsonl", "numeric_consistency.json")] + [
        PAIRED / name for name in ("run_config.json", "new_distributions.jsonl", "source_input_audit.csv", "summary.json",
                                  "paired_comparison.csv", "numeric_consistency.json")] + [
        SOURCES / "candidate_affinity.csv", SOURCES / "literature_table2.csv", SOURCES / "sources/PMC8463443.html",
        SCREEN / "skempi_index_screen.json", SCREEN / "raw/skempi_v2.csv",
        ROOT / "output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json"]
    config = {"audit_id": "EXISTING_CONTEXT_DISAGREEMENTS_20261022",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "authorization": "User 行，继续尝试吧 authorizes existing-score disagreement mining and local record matching only",
        "new_model_inputs": 0, "new_structure_jobs": 0, "new_reward_calls": 0,
        "downloads": 0, "training": 0, "numeric_tolerance_nats": 1e-4,
        "rules": "Reuse existing d sign and two-control H consensus; distinguish strict opposition from H abstention",
        "derived_arithmetic": "CPU d/H recomputation from saved distributions; no new model calls",
        "matching": "Exact full parent AND full variant AND receptor; no truncation or L10A-background transfer",
        "local_source_scope": ["Li2010 Table1/3 prior row audit", "Li2021 MDM2 L10A forward pairs",
                               "SKEMPI named 3EQS_A_B / 3LNZ_A_B / 3EQY_A_C only"],
        "posthoc_development_only": True, "primary_formal_metrics": None, "cpu_wall_seconds_max": 120.,
        "disk_bytes_max": 16777216, "preservation_before": preserved,
        "input_files": [{"path": str(p), "sha256": sha256(p)} for p in input_paths],
        "code_files": [{"path": str(p), "sha256": sha256(p)} for p in code_paths()]}
    write_json(HERE / "audit_config.json", config)
    (HERE / "audit_config.sha256").write_text(sha256(HERE / "audit_config.json") + "  audit_config.json\n")
    print({"status": "CPU_AUDIT_FROZEN", "new_model_inputs": 0, "inputs": len(input_paths),
           "config_sha256": sha256(HERE / "audit_config.json")}, flush=True)


def validate() -> dict:
    config = read_json(HERE / "audit_config.json")
    if sha256(HERE / "audit_config.json") != (HERE / "audit_config.sha256").read_text().split()[0]:
        raise ValueError("Audit config changed")
    for item in config["input_files"] + config["code_files"]:
        if sha256(Path(item["path"])) != item["sha256"]:
            raise ValueError(f"Frozen file changed: {item['path']}")
    return config


def experimental_inventory(target: str) -> list[dict]:
    experiments = []
    for row in read_csv(PAIRED / "source_input_audit.csv"):
        parent, variant = row["parent_sequence"], row["variant_sequence"]
        if len(parent) != len(variant) or len(parent) != 12 or parent == variant:
            continue
        changed = [i for i, (wt, alt) in enumerate(zip(parent, variant), 1) if wt != alt]
        if len(changed) != 1:
            continue
        experiments.append({"experiment_id": "Li2010_" + row["mutation_id"], "source": "Li2010",
            "source_path": str(PAIRED / "source_input_audit.csv"), "parent_sequence": parent, "variant_sequence": variant,
            "site": changed[0], "wt": parent[changed[0]-1], "alt": variant[changed[0]-1],
            "receptor_sequence": target, "receptor_construct": "MDM2(25-109)", "assay": row["assay"],
            "point_ratio": row["ratio_from_means"], "point_direction": row["point_direction"],
            "source_decision": row["decision"], "termini": row["termini"], "formal_label": "HOLD_NOT_RELEASED",
            "mutation_id": row["mutation_id"]})
    candidates = read_csv(SOURCES / "candidate_affinity.csv")
    parent_row = next(row for row in candidates if row["peptide_id"] == "MDM2_L10A")
    published = read_csv(SOURCES / "literature_table2.csv")
    parent_kd = float(published[0]["mdm2_kd_reported_nM"].split("±")[0].strip())
    if published[0]["sequence"] != PMI_L10A or parent_kd != float(parent_row["kd_nM"]):
        raise ValueError("Published Table2 parent differs from curated metadata")
    for row in published[1:]:
        parent, variant = PMI_L10A, row["sequence"]
        changes = [i for i, (wt, alt) in enumerate(zip(parent, variant), 1) if wt != alt]
        if len(parent) != 12 or len(variant) != 12 or len(changes) != 1:
            raise ValueError("Curated L10A forward panel no longer consists of single substitutions")
        kd = float(row["mdm2_kd_reported_nM"].split("±")[0].strip())
        ratio = kd / parent_kd
        mutation = f"{parent[changes[0]-1]}{changes[0]}{variant[changes[0]-1]}"
        experiments.append({"experiment_id": "Li2021_L10A_" + mutation, "source": "Li2021",
            "source_path": str(SOURCES / "literature_table2.csv"), "parent_sequence": parent, "variant_sequence": variant,
            "site": changes[0], "wt": parent[changes[0]-1], "alt": variant[changes[0]-1],
            "receptor_sequence": target, "receptor_construct": parent_row["receptor_construct_reported"], "assay": parent_row["assay"],
            "point_ratio": ratio, "point_direction": "gain" if ratio <= .5 else "loss" if ratio >= 2 else "ambiguous",
            "source_decision": "HOLD_L10A_BACKGROUND_ERROR_SEMANTICS_AND_TERMINI",
            "termini": parent_row["N_terminal"] + "; " + parent_row["C_terminal"], "formal_label": "HOLD_NOT_RELEASED",
            "mutation_id": "L10A_background_" + mutation})
    return experiments


def named_skempi_lookup() -> list[dict]:
    """Look up only the three previously identified groups; never rescreen leads."""
    names = {"3EQS_A_B", "3LNZ_A_B", "3EQY_A_C"}
    rows = [r for r in read_csv(SCREEN / "raw/skempi_v2.csv", ";") if r["#Pdb"] in names]
    result = []
    for row in rows:
        code = row["Mutation(s)_PDB"]
        match = re.fullmatch(r"([A-Z])([A-Z])(\d+)([A-Z])", code)
        if not match or row["Reference"] != "20226197":
            raise ValueError("Named SKEMPI lookup has an unexpected mutation/source")
        wt, chain, position, alt = match.groups()
        site = int(position)
        parent = "TSFAEYWALLSP" if row["#Pdb"] == "3LNZ_A_B" else PMI
        if parent[site-1] != wt:
            raise ValueError("SKEMPI mutation does not match its source parent")
        result.append({"pdb_chain_pair": row["#Pdb"], "mutation": code, "parent_sequence": parent,
            "variant_sequence": variant_sequence(parent, site, alt), "reference": row["Reference"],
            "method": row["Method"], "kd_parent_M": row["Affinity_wt (M)"], "kd_variant_M": row["Affinity_mut (M)"],
            "status": "EXCLUDE_MDMX_RECEPTOR" if row["#Pdb"] == "3EQY_A_C" else
                      "REVERSE_EDGE_SAME_LI2010_N8A_PAIR_NOT_NEW_SAMPLE" if row["#Pdb"] == "3LNZ_A_B" else
                      "DUPLICATE_LI2010_FORWARD_NOT_NEW_LABEL"})
    if len(result) != 22:
        raise ValueError("Named SKEMPI row count changed")
    return result


def worker() -> None:
    config = validate()
    old = read_json(PROBE / "run_config.json")
    paired = read_json(PAIRED / "run_config.json")
    for directory in (PROBE, PAIRED):
        if not read_json(directory / "numeric_consistency.json")["passed"]:
            raise ValueError("Existing distributions failed numerical consistency")
    weight = next(a["sha256"] for a in old["assets"]["assets"] if a["path"].endswith("pytorch_model.bin"))
    if weight != paired["weight_sha256"] or old["contexts"] != paired["contexts"]:
        raise ValueError("Cannot merge different checkpoints or receptor contexts")
    if ((old["precision"], old["attn_implementation"], old["tf32"], old["numeric_tolerance_nats"]) !=
            ("FP32", "eager", False, 1e-4) or paired["precision_key"] != "FP32_eager_TF32off" or
            paired["numeric_tolerance_nats"] != 1e-4):
        raise ValueError("Stored precision/numerical definitions differ")
    contexts = {r["id"]: r["sequence"] for r in old["contexts"]}
    target = old["contexts"][0]["sequence"]
    target_id = old["contexts"][0]["id"]
    old_parents = {p["parent_id"]: p["full_sequence"] for p in old["inputs"]}
    base = {}
    lineage = []
    for file_path, config_id in ((PROBE / "raw_masked_distributions.jsonl", "probe"),
                                 (PAIRED / "new_distributions.jsonl", "paired")):
        for record in jsonl(file_path):
            if record["phase"] != "base":
                continue
            parent = old_parents[record["parent_id"]] if config_id == "probe" else paired["parent_sequence"]
            key = distribution_key(weight, contexts[record["context_id"]], parent, record["site"])
            if key in base:
                raise ValueError("Duplicate base distribution rather than a documented reuse")
            if record["wt"] != parent[record["site"]-1] or any(not math.isfinite(v) for v in record["canonical_log_probs"].values()):
                raise ValueError("Invalid stored distribution")
            base[key] = record
            lineage.append({"source_file": str(file_path), "source_job_id": record["job_id"],
                            "parent_sequence": parent, "site": record["site"], "context_id": record["context_id"]})
    write_csv("distribution_lineage.csv", lineage)
    experiments = experimental_inventory(target)
    write_csv("experimental_pair_inventory.csv", experiments)
    skempi = named_skempi_lookup()
    write_csv("skempi_named_lookup.csv", skempi)
    sites = sorted({(key[2], key[3]) for key in base})
    if len(sites) != 19 or len(base) != 57:
        raise ValueError("Frozen existing-input scope differs from 19 sites x 3 contexts")
    all_predictions = []
    exact_matches = []
    near_matches = []
    original_pool = {p["full_sequence"]: set(p["eligible_positions"]) for p in old["inputs"]}
    history_exclusions = {("ETFSDLWKLLPEN", 6, "A"), (PMI, 8, "A")}
    for parent, site in sites:
        wt = parent[site-1]
        probs = [base[distribution_key(weight, c["sequence"], parent, site)]["canonical_log_probs"] for c in old["contexts"]]
        for alt in AA:
            if alt == wt:
                continue
            values = contrast(*probs, wt, alt)
            dp = predict_d(values["d_target"])
            hp = predict_h(values["h_control_1"], values["h_control_2"])
            parent_name = "PMI12" if parent == PMI else "P53_13MER"
            row = {"prediction_id": f"{parent_name}_{wt}{site}{alt}", "parent_sequence": parent,
                   "variant_sequence": variant_sequence(parent, site, alt), "site": site, "wt": wt, "alt": alt,
                   "receptor_sequence": target, **values, "d_prediction": dp, "h_prediction": hp,
                   "disagreement_type": disagreement_type(dp, hp), "raw_H_prediction": predict_d(values["h"]),
                   "original_283_generation_pool": site in original_pool.get(parent, set()) and (parent, site, alt) not in history_exclusions,
                   "historical_promotion_excluded": (parent, site, alt) in history_exclusions or (parent == PMI and site == 7 and alt == "A"),
                   "new_design_recommendation": False, "exact_experimental_matches": 0,
                   "eligible_conditional_development_matches": 0, "formal_labels": 0}
            for experiment in experiments:
                if exact_pair(row, experiment):
                    row["exact_experimental_matches"] += 1
                    usable = experiment["source_decision"] == "CONDITIONAL_DEVELOPMENT_ONLY"
                    row["eligible_conditional_development_matches"] += int(usable)
                    exact_matches.append({"prediction_id": row["prediction_id"], "experiment_id": experiment["experiment_id"],
                        "parent_sequence": parent, "variant_sequence": row["variant_sequence"], "mutation_id": experiment["mutation_id"],
                        "disagreement_type": row["disagreement_type"], "d_prediction": dp, "h_prediction": hp,
                        "point_ratio": experiment["point_ratio"], "point_direction": experiment["point_direction"],
                        "source_decision": experiment["source_decision"], "conditional_development_eligible": usable,
                        "formal_label": "HOLD_NOT_RELEASED"})
                else:
                    reason = near_pair_reason(row, experiment)
                    if reason:
                        near_matches.append({"prediction_id": row["prediction_id"], "experiment_id": experiment["experiment_id"],
                            "model_parent": parent, "experimental_parent": experiment["parent_sequence"],
                            "model_variant": row["variant_sequence"], "experimental_variant": experiment["variant_sequence"],
                            "disagreement_type": row["disagreement_type"], "reason": reason,
                            "label_transferred": False})
            all_predictions.append(row)
    disagreements = [r for r in all_predictions if r["disagreement_type"] in DISAGREEMENT_TYPES]
    write_csv("all_predictions.csv", all_predictions)
    write_csv("disagreements.csv", disagreements, list(all_predictions[0]))
    write_csv("exact_experimental_matches.csv", exact_matches, ["prediction_id", "experiment_id"])
    write_csv("near_matches_hold.csv", near_matches, ["prediction_id", "experiment_id", "reason"])
    source_scope = [
        {"source": "Li2010 Table1/3", "local_pair_records": sum(r["source"] == "Li2010" for r in experiments),
         "role": "Original source pairs; no formal labels released"},
        {"source": "Li2021 L10A MDM2 forward", "local_pair_records": sum(r["source"] == "Li2021" for r in experiments),
         "role": "L10A parent mismatch; undefined error semantics and termini remain HOLD"},
        {"source": "SKEMPI 3EQS/3LNZ", "local_pair_records": sum(r["pdb_chain_pair"] != "3EQY_A_C" for r in skempi),
         "role": "Same Li2010 source; forward duplicates and reverse edge are not new samples"},
        {"source": "SKEMPI 3EQY", "local_pair_records": sum(r["pdb_chain_pair"] == "3EQY_A_C" for r in skempi),
         "role": "MDMX receptor excluded"}]
    write_csv("source_scope_audit.csv", source_scope)
    by_parent = {}
    for parent in sorted({r["parent_sequence"] for r in all_predictions}):
        group = [r for r in all_predictions if r["parent_sequence"] == parent]
        by_parent[parent] = {"n": len(group), "types": dict(Counter(r["disagreement_type"] for r in group)),
                             "original_generation_pool": sum(r["original_283_generation_pool"] for r in group)}
    exact_disagreements = [r for r in exact_matches if r["disagreement_type"] in DISAGREEMENT_TYPES]
    eligible_disagreements = [r for r in exact_disagreements if r["conditional_development_eligible"]]
    summary = {"status": "EXISTING_SCORE_DISAGREEMENT_AUDIT_COMPLETED", "new_model_inputs": 0,
        "new_structure_jobs": 0, "new_reward_calls": 0, "downloads": 0, "training": 0,
        "unique_base_distributions": len(base), "unique_mask_sites": len(sites), "non_wt_substitutions": len(all_predictions),
        "original_generation_pool": sum(r["original_283_generation_pool"] for r in all_predictions),
        "disagreement_counts": dict(Counter(r["disagreement_type"] for r in all_predictions)),
        "by_parent": by_parent, "disagreement_subset_size": len(disagreements),
        "exact_experimental_matches_all": len(exact_matches), "exact_experimental_matches_in_disagreements": len(exact_disagreements),
        "eligible_conditional_development_matches_in_disagreements": len(eligible_disagreements),
        "formal_eligible_matches_in_disagreements": 0, "primary_formal_metrics": None,
        "near_match_holds": dict(Counter(r["reason"] for r in near_matches)),
        "skempi_named_rows": len(skempi), "skempi_status_counts": dict(Counter(r["status"] for r in skempi)),
        "source_scope": source_scope, "posthoc_development_only": True,
        "selection_warning": "Score-selected disagreement subset is posthoc development; never a blind or representative benchmark",
        "decision": "HOLD_NO_ADJUDICATED_DISAGREEMENT" if not eligible_disagreements else "CONDITIONAL_DEVELOPMENT_ONLY_REVIEW_MATCHES",
        "affinity_validation": "HOLD", "closed_loop": "HOLD"}
    write_json(HERE / "summary.json", summary)
    no_model_modules = not any(name in sys.modules for name in ("torch", "transformers", "jax", "esm"))
    write_json(HERE / "runtime_environment.json", {"python": sys.executable, "version": sys.version,
               "model_modules_not_imported": no_model_modules, "new_model_inputs": 0})
    if not no_model_modules:
        raise ValueError("A model runtime module was unexpectedly imported during CPU-only audit")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


def run() -> None:
    config = validate()
    if (HERE / "audit_worker.log").exists():
        raise FileExistsError("An audit attempt already exists")
    result = run_bounded([sys.executable, "-B", str(Path(__file__).resolve()), "worker"], HERE,
                          config["cpu_wall_seconds_max"], config["disk_bytes_max"], HERE / "audit_worker.log")
    write_json(HERE / "budget_result.json", result)
    print(json.dumps(result), flush=True)
    if result["status"] != "completed":
        raise RuntimeError("CPU audit stopped; preserve partial outputs")


def deliver() -> None:
    config = validate()
    summary = read_json(HERE / "summary.json")
    budget = read_json(HERE / "budget_result.json")
    tests = subprocess.run([sys.executable, "-B", str(HERE / "test_audit.py")], capture_output=True, text=True)
    (HERE / "test_after_run.log").write_text(tests.stdout + tests.stderr, encoding="utf-8")
    after = history()
    all_rows = read_csv(HERE / "all_predictions.csv")
    disagreements = read_csv(HERE / "disagreements.csv")
    exact = read_csv(HERE / "exact_experimental_matches.csv")
    exact_disagreements = [r for r in exact if r["disagreement_type"] in DISAGREEMENT_TYPES]
    near = read_csv(HERE / "near_matches_hold.csv")
    counter = summary["disagreement_counts"]
    if not disagreements:
        interpretation = "本次已有分布未发现按冻结规则可区分的方向/拒绝变化，不能据此评价上下文修正增益。"
    elif summary["eligible_conditional_development_matches_in_disagreements"] == 0:
        interpretation = "已有读数确有不同方向或拒绝，但本地记录中尚无满足条件开发定义的分歧标签，不能判定谁更可靠。"
    else:
        interpretation = "分歧中有条件开发配对，必须按完整匹配和误差限制逐条核对；这些仍不构成正式或独立验证。"
    checks = [
        {"name": "model_free_tests", "pass": tests.returncode == 0},
        {"name": "frozen_inputs_and_code", "pass": True},
        {"name": "prior_manifests_unchanged", "pass": after == config["preservation_before"] and not any(p["mismatches"] for p in after)},
        {"name": "361_rows_19_sites_no_duplicate_mutants", "pass": len(all_rows) == 361 and len({(r["parent_sequence"], r["site"]) for r in all_rows}) == 19 and
         len({(r["parent_sequence"], r["variant_sequence"]) for r in all_rows}) == 361},
        {"name": "original_283_pool_preserved", "pass": sum(r["original_283_generation_pool"] == "True" for r in all_rows) == 283},
        {"name": "all_prediction_types_retained", "pass": sum(counter.values()) == 361 and len(disagreements) == summary["disagreement_subset_size"]},
        {"name": "no_near_label_transfer", "pass": all(r["label_transferred"] == "False" for r in near)},
        {"name": "formal_metrics_not_released", "pass": summary["primary_formal_metrics"] is None and summary["formal_eligible_matches_in_disagreements"] == 0},
        {"name": "p27_not_in_prediction_sites", "pass": not any(r["parent_sequence"] == "ETFSDLWKLLPEN" and r["site"] == "11" for r in all_rows)},
        {"name": "no_model_or_structure_or_reward_calls", "pass": all(summary[k] == 0 for k in ("new_model_inputs", "new_structure_jobs", "new_reward_calls", "downloads", "training"))},
        {"name": "no_model_runtime_imports", "pass": read_json(HERE / "runtime_environment.json")["model_modules_not_imported"]},
        {"name": "bounded_cpu_completion", "pass": budget["status"] == "completed" and budget["wall_seconds"] < 120 and directory_bytes(HERE) < 16777216},
    ]
    lines = ["# 已有H/d分歧与实验配对审计", "",
        f"结论：{summary['decision']}。本轮只处理已有读数，新增模型前向、结构预测、GraphPep/FoldX/ESM2调用、下载和训练均为0；H/d由已有分布作CPU算术复算。", "",
        "## 实际分歧", "",
        f"两轮已完成推理共有57条唯一基础分布，覆盖19个位点，读取每个位点19个非WT替换后共{len(all_rows)}项。",
        "p53输入是13mer，PMI是12mer；位点来源/完整序列/受体/checkpoint均保留。原283个生成候选池单独标记；",
        "后来仅用于评价的锚位及历史排除没有变成设计推荐。CSV中完整变体序列只用于实验身份匹配，不代表执行了新生成或推荐。", "",
        "H/d沿用先前固定方向规则和1e-4 nats容差，未根据本轮标签调整。H的两控制不同向时为拒绝，不称方向反转。", "",
        "|模型亲本|已有替换项|同向|严格反向|H拒绝|d平局而H有方向|", "|---|---:|---:|---:|---:|---:|"]
    for parent, item in summary["by_parent"].items():
        counts = item["types"]
        lines.append(f"|{parent}|{item['n']}|{counts.get('same_direction',0)}|{counts.get('strict_opposite',0)}|"
                     f"{counts.get('H_abstains',0)}|{counts.get('d_tie_H_directional',0)}|")
    lines.extend(["", f"合计分歧子集{len(disagreements)}项。完整361项和所有拒绝均保留，不报告挑选分歧子集后的总体准确率。", "",
        "## 本地实验匹配", "",
        f"完整亲本与完整变体均精确匹配的实验配对有{summary['exact_experimental_matches_all']}项；分歧子集中有{len(exact_disagreements)}项。",
        f"其中满足上轮条件开发进入定义的分歧有{summary['eligible_conditional_development_matches_in_disagreements']}项；正式合格分歧仍为0，正式主指标=null。", "",
        "|分歧预测|原实验|点Kd比值|d方向|H方向|来源判定|", "|---|---|---:|---|---|---|"])
    for row in exact_disagreements:
        lines.append(f"|{row['prediction_id']}|{row['mutation_id']}|{float(row['point_ratio']):.6g}|"
                     f"{row['d_prediction']}|{row['h_prediction']}|{row['source_decision']}|")
    if not exact_disagreements:
        lines.append("|无完整配对实验命中|—|—|—|—|HOLD|")
    lines.extend(["", "点Kd变化未达冻结两倍门的项目不能被升级为增强/减弱标签。近中性读数既不证明H纠错，也不证明H引错。",
        "此前七对PMI条件开发例的H/d仍全部同向；不能从这些7/7结果推导上下文修正有额外收益。", "",
        "## 不允许的标签移植及重复证据", "",
        f"近似身份匹配HOLD：{json.dumps(summary['near_match_holds'], ensure_ascii=False)}。",
        "p53 13mer不能裁成12mer继承实验Kd；PMI-L10A与PMI-WT是不同背景，不能仅按N8A等突变名拼接。",
        "即使最终variant偶然相同，也不允许把另一个parent的实验差值移过来。2021表的端基/误差语义缺项保持原判定。",
        f"只查既有SKEMPI索引指定的22行：{json.dumps(summary['skempi_status_counts'], ensure_ascii=False)}。",
        "3EQS与3LNZ条目仍是同一Li2010来源；N8A反向边不是新的减弱样本。3EQY是MDMX，不能移作MDM2证据。",
        "此处只核对已知条目，没有重新筛全库、下载附件或将重复数据包装成新进展。", "",
        "## 解释与有限下一步", "",
        interpretation,
        "该分歧清单按已经看过的分数筛出，属于后验开发诊断；后续补到标签仍不能把它称盲测或独立测试。",
        "本轮不扩大H，不重新计算相同PMI面板，不为未知标签启动折叠。保留固定分歧清单及身份缺项，下一次若有精确配对实验记录，",
        "可按原规则核对H纠错或引错；在此之前没有依据选择主reward。PDB能提供结构与映射，不能填补这些配对的实验Kd。",
        "正式亲和力验证和设计闭环继续HOLD，不改变旧结构质量门、0.05平局带或N8A不确定性限制。", "",
        "## 可复算产物", "",
        "all_predictions.csv保留361项；disagreements.csv保留反向与拒绝变化；exact_experimental_matches.csv记录完整配对命中；",
        "near_matches_hold.csv保留输入不匹配的近似命中；experimental_pair_inventory.csv与skempi_named_lookup.csv保留有限来源；",
        "distribution_lineage.csv连接原始job；audit_config.json/sha256、预算、测试与verification.json记录本轮执行。"])
    (HERE / "report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(HERE / "verification.json", {"checks": checks, "failed_count": sum(not c["pass"] for c in checks),
        "preservation": after, "scientific_scope": "Posthoc existing-score/local-record audit; no new affinity validation"})
    if any(not c["pass"] for c in checks):
        raise ValueError("Audit delivery verification failed")
    files = [{"path": str(p.relative_to(HERE)), "bytes": p.stat().st_size, "sha256": sha256(p)}
             for p in sorted(HERE.rglob("*")) if p.is_file() and p.name not in
             ("artifact_manifest.json", "artifact_manifest.sha256")]
    write_json(HERE / "artifact_manifest.json", {"scope": "All files except this manifest and detached checksum", "files": files})
    (HERE / "artifact_manifest.sha256").write_text(sha256(HERE / "artifact_manifest.json") + "  artifact_manifest.json\n")
    print({"checks": len(checks), "failed": sum(not c["pass"] for c in checks), "artifact_files": len(files),
           "output_bytes": directory_bytes(HERE)}, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "run", "worker", "deliver"))
    globals()[parser.parse_args().mode]()
