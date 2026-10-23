"""Post-run checks and Chinese development report; no model calls."""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import json
import subprocess
import sys

from run_development import HERE, OLD, TOL, preserve, validate
from probe_core import directory_bytes, read_json, sha256, write_json


def main() -> None:
    config = validate()
    budget = read_json(HERE / "budget_result.json")
    tests = subprocess.run([sys.executable, "-B", str(HERE / "test_paired_core.py")], capture_output=True, text=True)
    (HERE / "test_after_run.log").write_text(tests.stdout + tests.stderr, encoding="utf-8")
    history = preserve()
    checks = [
        {"name": "model_free_tests", "pass": tests.returncode == 0},
        {"name": "frozen_input_code_model_hashes", "pass": True},
        {"name": "historical_artifacts_preserved", "pass": history == config["preservation_before"] and
         not any(item["mismatches"] for item in history)},
        {"name": "wall_budget", "pass": budget["wall_seconds"] <= 120.},
        {"name": "disk_budget", "pass": directory_bytes(HERE) < 67108864},
    ]
    if not (HERE / "summary.json").exists():
        progress = read_json(HERE / "progress.json") if (HERE / "progress.json").exists() else {}
        lines = ["# PMI开发对照：停止记录", "",
                 f"本轮worker状态：{budget['status']}，GPU进程墙钟{budget['wall_seconds']:.3f}秒。", "",
                 f"进度：{progress}。原始失败、局部读数与日志保留；没有重试或放宽容差。", "",
                 "本轮不报告方向比较或算法增益。正式v2标签、亲和力验证和闭环均继续HOLD。"]
    else:
        summary = read_json(HERE / "summary.json")
        runtime = read_json(HERE / "runtime.json")
        numeric = summary["numeric"]
        with (HERE / "paired_comparison.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        with (HERE / "source_input_audit.csv").open(newline="") as handle:
            source_audit = list(csv.DictReader(handle))
        reused = [json.loads(line) for line in (HERE / "reused_distributions.jsonl").read_text().splitlines()]
        new = [json.loads(line) for line in (HERE / "new_distributions.jsonl").read_text().splitlines()]
        all_base = {(r["site"], r["context_id"]): r for r in reused + new if r["phase"] == "base"}
        context_ids = [c["id"] for c in config["contexts"]]
        pairs = {p["mutation_id"]: p for p in config["pairs"]}
        recalculated = True
        predictions = True
        mapping = True
        for row in rows:
            pair = pairs[row["mutation_id"]]
            site = int(row["site"])
            wt, alt = pair["parent_sequence"][site-1], pair["variant_sequence"][site-1]
            odds = []
            for cid in context_ids:
                logp = all_base[(site, cid)]["canonical_log_probs"]
                odds.append(logp[alt] - logp[wt])
            h1, h2 = odds[0]-odds[1], odds[0]-odds[2]
            calculated = {"d_target": odds[0], "d_control_1": odds[1], "d_control_2": odds[2],
                          "h_control_1": h1, "h_control_2": h2, "h": odds[0] - sum(odds[1:])/2}
            recalculated &= all(abs(float(row[key])-value) < 1e-12 for key, value in calculated.items())
            expected_d = "gain" if odds[0] > TOL else "loss" if odds[0] < -TOL else "tie"
            expected_h = "gain" if h1 > TOL and h2 > TOL else "loss" if h1 < -TOL and h2 < -TOL else "abstain"
            predictions &= row["d_prediction"] == expected_d and row["h_prediction"] == expected_h
            mapping &= (site, row["wt"], row["alt"], row["point_direction"]) == (pair["site"], wt, alt, pair["point_direction"])
        old_records = {r["job_id"]: r for r in
                       [json.loads(line) for line in (OLD / "raw_masked_distributions.jsonl").read_text().splitlines()]}
        cache_values_exact = all(r["canonical_log_probs"] == old_records[r["reused_from_job_id"]]["canonical_log_probs"]
                                 and r["full_cache_key"][2] == config["parent_sequence"] for r in reused)
        comparison = summary["comparison"]
        checks.extend([
            {"name": "new_19_reused_9_base_21", "pass": len(new) == 19 and len(reused) == 9 and len(all_base) == 21 and
             runtime["logical_inputs_completed"] == 19 and budget["status"] == "completed"},
            {"name": "seven_repeat_numeric_gate", "pass": numeric["passed"] and len(numeric["checks"]) == 7 and
             max(numeric["max_logp_delta"], numeric["max_d_delta"]) <= TOL},
            {"name": "all_seven_pairs_no_denominator_drops", "pass": len(rows) == 7 and len({r["mutation_id"] for r in rows}) == 7 and
             comparison["n"] == 7 and sum(comparison[k] for k in ("both_correct", "h_only_correct", "d_only_correct", "neither_correct")) == 7},
            {"name": "direct_raw_distribution_score_recalculation", "pass": recalculated},
            {"name": "independent_prediction_rule_check", "pass": predictions},
            {"name": "sequence_site_and_source_mapping", "pass": mapping},
            {"name": "reused_distribution_values_exact", "pass": cache_values_exact},
            {"name": "formal_metrics_and_gain_ci_not_released", "pass": summary["primary_formal_metrics"] is None and
             summary["formal_method_difference_95_CI"] is None and summary["formal_labels_released"] == 0},
            {"name": "p53_and_p27_not_in_jobs", "pass": all(j["parent_id"] == "PMI_WT_12MER" and
             j["mutation_id"].startswith("PMI_") for j in config["jobs"])},
            {"name": "fp32_eval_tf32_off", "pass": runtime["parameter_dtype"] == "float32" and runtime["eval"] and
             not runtime["tf32_matmul"] and not runtime["tf32_cudnn"]},
        ])
        labels = {"gain": "增强", "loss": "减弱", "tie": "平局", "abstain": "拒绝", "failure": "失败"}
        point_gains = comparison["by_direction"]["gain"]
        point_losses = comparison["by_direction"]["loss"]
        n8 = next(row for row in rows if row["mutation_id"] == "PMI_N8A")
        counts = Counter(row["decision"] for row in source_audit)
        lines = ["# PMI条件开发对照：上下文修正H与原始d", "",
            f"结论：{summary['decision']}。本轮按已冻结规则完成，d符合{comparison['d_correct']}/7，",
            f"H（控制一致方向）符合{comparison['h_correct']}/7；H拒绝/失败{comparison['h_abstain_or_failure']}/7，均保留分母。",
            "这些是对文献点效应方向的开发核对，不是正式方向准确率或独立验证；正式主指标和正式方法增益95%区间均为null。", "",
            "## 实际执行与输入对应", "",
            "用户明确“OK，开始吧”启动此前准备的七对PMI方案。只读本地Li 2010原文Table 1/3，未重新下载或扩展来源。",
            "完整亲本TSFAEYWNLLSP（12mer），受体MDM2(25–109)85aa，同一SPR competition assay。",
            "原文明确肽未N-acetylated、未C-amidated，模型中的标准自由N/C端与该端基身份对应；模型概率不复现溶液条件或直接测量Kd。",
            "七个实验替换位点可含锚点，因为本任务是评价已知突变，未修改旧候选生成的锚点保护规则。",
            f"21条基础分布中精确复用9条，新计算12条基础加7条目标上下文逆序重复，共{len(new)}条新逻辑输入。",
            "复用键含完整受体序列、完整肽序列、位点、checkpoint、精度和固定评分定义，不使用peptide-only PPL缓存。",
            f"FP32/eval/no gradients/TF32 off/eager，batch≤2；GPU进程墙钟{budget['wall_seconds']:.3f}/120秒（含加载与进程内核验），",
            f"峰值CUDA allocated {runtime['peak_cuda_allocated_bytes']/2**30:.3f} GiB，reserved {runtime['peak_cuda_reserved_bytes']/2**30:.3f} GiB。",
            f"七次重复最大logp差{numeric['max_logp_delta']:.9g}、最大d差{numeric['max_d_delta']:.9g} nats，均≤1e-4。",
            "没有训练、采样完整子代、预测新结构或调用GraphPep/FoldX/ESM2；没有重跑P27A。", "",
            "## 冻结方法", "",
            "d=logp(替换残基|目标受体,单点mask亲本)−logp(原残基|相同输入)。d>1e-4预测增强，d<−1e-4预测减弱，否则tie。",
            "H=d_target−(d_shuffle1+d_shuffle2)/2；两控制修正分别为H1=d_target−d_shuffle1、H2=d_target−d_shuffle2。",
            "H1/H2均>1e-4预测增强，均<−1e-4预测减弱，否则拒绝；拒绝保留在七对分母。两个固定打乱控制不是实验非结合受体。",
            "1e-4只是nats数值容差，无生物学效应量含义；GraphPep/FoldX既有0.05平局带保持不变。",
            "原来的增强建议filter仅作为诊断列，不用于筛测试集。raw H符号固定为辅助消融，不按结果选方法。", "",
            "## 逐对结果（有条件开发核对）", "",
            "|突变|Kd比值（均值重算）|文献点效应|d(nats)|d方向|H(nats)|H方向|",
            "|---|---:|---|---:|---|---:|---|"]
        for row in rows:
            lines.append(f"|{row['mutation_id']}|{float(row['ratio_from_means']):.6g}|{labels[row['point_direction']]}|"
                         f"{float(row['d_target']):.6f}|{labels[row['d_prediction']]}|{float(row['h']):.6f}|{labels[row['h_prediction']]}|")
        lines.extend(["", "|点效应类别|对数|d符合数|H符合数|H拒绝/失败数|",
            "|---|---:|---:|---:|---:|",
            f"|增强|{point_gains['n']}|{point_gains['d_correct']}|{point_gains['h_correct']}|{point_gains['h_abstain_or_failure']}|",
            f"|减弱|{point_losses['n']}|{point_losses['d_correct']}|{point_losses['h_correct']}|{point_losses['h_abstain_or_failure']}|", "",
            f"成对核对：两者均符合{comparison['both_correct']}；仅H符合{comparison['h_only_correct']}；",
            f"仅d符合{comparison['d_only_correct']}；两者均不符合{comparison['neither_correct']}。",
            f"raw H只按符号的固定辅助消融为{summary['raw_H_sign_diagnostic_correct']}/7；不取代主表的控制一致规则。", "",
            "## 标签与不确定性限制", "",
            "七对只有一个亲本、一个增强点效应，且家族/旧反例已暴露。PMI正文报告±SD与至少三次独立测量，逐肽精确n、原始重复和协方差未恢复。",
            "敏感性列假设n=3，使用t(df=2)和log-ratio delta method，分别组合独立误差与相对SD之和的最坏相关情景。",
            "这些只是有假设的约95%敏感性范围，不是来源CI或保证保守的置信区间；N8A的变体CV约43%，线性近似也有限。",
            f"N8A独立误差情景比值范围[{float(n8['independent_sensitivity_lower']):.6g}, {float(n8['independent_sensitivity_upper']):.6g}]；",
            f"最坏相关情景[{float(n8['worst_correlation_sensitivity_lower']):.6g}, {float(n8['worst_correlation_sensitivity_upper']):.6g}]跨1，不能称稳健正式增强标签。",
            "因此全部正式v2标签仍未释放，正式主指标=null。不计算伪精确的正式成对增益95%区间；逐对差异和类别结果完整报告。", "",
            "## 来源审计与HOLD", "",
            f"原两表共保存{len(source_audit)}个HTML网格行，其中{counts['NONDATA_SEPARATOR']}个为空分隔行；所有决定见source_input_audit.csv。",
            f"PMI七对为条件开发任务；p53满足两倍点效应的{counts['HOLD_REPLICATE_COUNT_NOT_ATTRIBUTABLE_TO_P53']}对因重复信息归属不足而HOLD，没有新推理。",
            "p53真正assay输入为12mer ETFSDLWKLLPE；旧13mer不能继承该标签。K24A=230/440≈0.523与E28A=240/440≈0.545不符合≤0.5，",
            "不能按印刷表四舍五入的0.5纳入增强。P27A排除；ND没有可证明方向的检测限，不伪造点Kd。", "",
            "## 决定与下一步", ""])
        if summary["development_paired_gain_observed"]:
            lines.append("本小面板出现H相对d的开发符合数增加，仅保留为后续验证线索；单一增强、单一背景及误差限制不足以宣称增益或选择主reward。")
        else:
            lines.append("本小面板未观察到H相对d的开发符合数增加，不支持扩大当前上下文修正版本。保留失败/相同结果，不改公式、控制或容差后再用同一组数据宣布成功。")
        lines.extend(["",
            "正式亲和力验证和新设计闭环继续HOLD。若以后提出新假设，应说明与本版的差别、可证伪对照和限定输入；已看过的面板仍只作开发。",
            "实验PDB仍可用于另行审计的给定pose结构评分；这不消除测量标签缺口，也不等于完整ΔG。本轮没有开启这条计算分支。", "",
            "## 可复算文件", "",
            "run_config.json/sha256冻结来源、输入、代码、方法及预算；new_distributions.jsonl保留19条新读数，reused_distributions.jsonl保留9条复用读数及来源job ID；",
            "paired_comparison.csv保留七对全部分数、方向、拒绝、条件敏感性与formal HOLD。summary.json、numeric_consistency.json、budget_result.json和verification.json为机器记录。",
            "旧proposal、六结构pilot、kernel ablation、首轮60输入probe的原始文件逐项核验，未改动；历史NOT_RUN/未勾选项按生成时点解释。"])
        write_json(HERE / "execution_record.json", {"completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "model_free_tests_passed": tests.returncode == 0, "new_inputs": len(new), "reused_inputs": len(reused),
            "pairs_retained": len(rows), "budget": budget, "formal_validation": "HOLD", "closed_loop": "HOLD"})
    (HERE / "report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(HERE / "verification.json", {"checks": checks, "failed_count": sum(not c["pass"] for c in checks),
        "preservation": history, "scope": "Technical and conditional-development verification only; no formal biological validation"})
    if any(not c["pass"] for c in checks):
        raise ValueError("Delivery verification failed")
    files = [{"path": str(p.relative_to(HERE)), "bytes": p.stat().st_size, "sha256": sha256(p)}
             for p in sorted(HERE.rglob("*")) if p.is_file() and p.name not in
             ("artifact_manifest.json", "artifact_manifest.sha256")]
    write_json(HERE / "artifact_manifest.json", {"scope": "All files except this manifest and detached checksum", "files": files})
    (HERE / "artifact_manifest.sha256").write_text(sha256(HERE / "artifact_manifest.json") + "  artifact_manifest.json\n")
    print({"checks": len(checks), "failed": sum(not c["pass"] for c in checks), "artifact_files": len(files),
           "output_bytes": directory_bytes(HERE), "manifest_sha256": sha256(HERE / "artifact_manifest.json")})


if __name__ == "__main__":
    main()
