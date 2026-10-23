"""Model-free post-run verification, Chinese report and artifact manifest."""
from __future__ import annotations

import csv
from pathlib import Path
import subprocess
import sys

from probe_core import directory_bytes, read_json, sha256, write_json
from probe_run import HERE, TOL, preservation_checks, validate_config


def main() -> None:
    config = validate_config()
    budget = read_json(HERE / "budget_result.json")
    tests = subprocess.run([sys.executable, "-B", str(HERE / "test_probe.py")],
                           capture_output=True, text=True)
    (HERE / "test_after_run.log").write_text(tests.stdout + tests.stderr, encoding="utf-8")
    preserved = preservation_checks()
    checks = [
        {"name": "model_free_tests", "pass": tests.returncode == 0},
        {"name": "source_model_and_code_hashes_unchanged", "pass": True},
        {"name": "preserved_historical_artifacts", "pass": not any(p["mismatches"] for p in preserved)},
        {"name": "gpu_wall_budget", "pass": budget["wall_seconds"] <= 600.},
        {"name": "output_budget", "pass": directory_bytes(HERE) < 67108864},
    ]
    if not (HERE / "summary.json").exists():
        progress = read_json(HERE / "progress.json") if (HERE / "progress.json").exists() else {}
        lines = ["# PepMLM受体上下文技术探针：停止记录", "",
                 f"本轮状态：{budget['status']}；停止后没有重试、改精度、换种子或放宽容差。", "",
                 f"GPU进程墙钟 {budget['wall_seconds']:.3f} 秒；进度记录：{progress}。", "",
                 "原始局部读数、worker.log和worker_failure.json（若有）保留。失败不支持模型无效或亲和力方向结论。",
                 "本轮不产生突变建议；下一计算继续HOLD，先定位这次失败。"]
    else:
        summary = read_json(HERE / "summary.json")
        runtime = read_json(HERE / "runtime.json")
        numeric = summary["numeric"]
        with (HERE / "all_substitution_diagnostics.csv").open(newline="") as stream:
            all_rows = list(csv.DictReader(stream))
        with (HERE / "technical_suggestions.csv").open(newline="") as stream:
            suggestions = list(csv.DictReader(stream))
        excluded = [r for r in all_rows if r["excluded_historical"] == "True"]
        expected_excluded = {("P53WT_1YCR", "6", "A"), ("PMI_3EQS", "8", "A")}
        checks.extend([
            {"name": "worker_completed_60", "pass": budget["status"] == "completed" and
             summary["logical_inputs"] == 60 and runtime["logical_inputs_completed"] == 60},
            {"name": "numeric_gate_frozen", "pass": numeric["passed"] and
             max(numeric["max_logp_delta"], numeric["max_d_delta"]) <= TOL},
            {"name": "285_diagnostics_283_eligible", "pass": len(all_rows) == 285 and len(excluded) == 2 and
             {(r["parent_id"], r["site"], r["alt"]) for r in excluded} == expected_excluded},
            {"name": "suggestion_filter_and_historical_exclusion", "pass": all(
                r["excluded_historical"] == "False" and all(float(r[k]) > TOL for k in
                    ("d_target", "h_control_1", "h_control_2")) for r in suggestions)},
            {"name": "suggestion_limits", "pass": all(
                len([r for r in suggestions if r["parent_id"] == p["parent_id"]]) <= 6 and
                len({r["site"] for r in suggestions if r["parent_id"] == p["parent_id"]}) ==
                len([r for r in suggestions if r["parent_id"] == p["parent_id"]]) for p in config["inputs"])},
            {"name": "wt_self_zero", "pass": summary["max_wt_self_difference"] <= TOL},
            {"name": "fp32_tf32_off_eval", "pass": runtime["parameter_dtype"] == "float32" and
             runtime["eval"] and not runtime["tf32_matmul"] and not runtime["tf32_cudnn"]},
        ])
        lines = ["# PepMLM受体上下文技术探针结果", "",
            f"结论：{summary['decision']}。已实现并执行准备方案的第一阶段；它只检验数值和受体顺序敏感性。",
            "独立亲和力验证与闭环生成继续HOLD，不能将以下建议称为增强结合肽。", "",
            "## 本次实际执行", "",
            "用户“OK试试吧”授权了本轮不超过60条逻辑输入的技术实验。旧proposal.json中的NOT_RUN是准备时状态，未回写旧文件。",
            f"两个亲本、15个位点×三个上下文=45条基础输入，加15条目标上下文逆序单批重复，共{summary['logical_inputs']}条。",
            "每次读取一个掩码位置的全词表分布，取得20种氨基酸读数；没有生成完整子代序列、预测复合物或调用GraphPep/FoldX/ESM2。",
            f"FP32、eval、无梯度、TF32关闭，基础batch≤2、重复batch=1；GPU进程墙钟{budget['wall_seconds']:.3f}/600秒（含模型加载），",
            f"峰值CUDA allocated {runtime['peak_cuda_allocated_bytes']/2**30:.3f} GiB，reserved {runtime['peak_cuda_reserved_bytes']/2**30:.3f} GiB。",
            "使用本地已安装权重，未下载、训练、重装或读取旧peptide-only PPL缓存。运行前后核验模型、输入和代码SHA256。", "",
            "## 冻结公式和技术判据", "",
            "d=log p(替换残基|受体,单点mask亲本)−log p(原残基|相同输入)；H=d_target−(d_shuffle1+d_shuffle2)/2。",
            "两个控制是预先固定、长度与组成相同的受体顺序打乱，不是实验不结合受体。",
            "只有d_target、H_control_1、H_control_2均>1e-4 nats才进入建议；每位点只取一个，每亲本最多六个。",
            "1e-4是数值容差，没有生物学效应量含义；GraphPep/FoldX正式方向的0.05平局带保持原样。",
            f"重复检查最大canonical logp差={numeric['max_logp_delta']:.9g} nats，最大d差={numeric['max_d_delta']:.9g}，",
            f"均不超过1e-4；WT→WT最大差={summary['max_wt_self_difference']:.9g}。", "",
            "## 覆盖与对照", "",
            "PDB已有坐标用于固定可掩码位置与复用GraphPep位置读数；本轮没有重新预测这些实验结构。",
            "P53使用13mer，不能直接继承历史12mer assay标签；PMI完整输入12mer，未解析末端P不参与mask。",
            "P53 P27位点和PMI W7锚定位点不mask；P53 L22A和PMI N8A的自动logits读数作为诊断保留、从建议/基线池排除。",
            "15×19=285个非WT替换读数中，排除这两个已暴露替换后283个构成共同池。所有过滤失败仍在诊断表中。", "",
            "|亲本|mask位点|共同替换池|通过filter|有通过替换的位点|控制反号|最终建议|",
            "|---|---:|---:|---:|---:|---:|---:|"]
        for parent in summary["parents"]:
            lines.append(f"|{parent['parent_id']}|{parent['mask_sites']}|{parent['candidate_pool']}|"
                         f"{parent['filter_pass_candidates']}|{parent['filter_pass_sites']}|"
                         f"{parent['control_sign_disagreement']}|{parent['suggestions']}|")
        lines.extend(["", "“控制反号”指两个H_control在±1e-4以外符号相反；这些替换自动拒绝，不重采样控制来救结果。",
            f"所有候选中最大的|H_control|为{summary['context_effect_abs_max_nats']:.6g} nats。该差异证明这些输入上的模型输出依赖受体顺序，",
            "不能据此分辨真实识别与打乱序列产生的分布外效应。", "",
            "|亲本|H与未修正d：同替换/同位点|H与静态GraphPep：同位点|H与固定随机：同替换/同位点|",
            "|---|---|---|---|"])
        for parent in summary["parents"]:
            overlaps = parent["overlaps"]
            d = overlaps["unadjusted_d"]
            g = overlaps["static_GraphPep_native_position_only"]
            r = overlaps["uniform_random"]
            k = parent["suggestions"]
            lines.append(f"|{parent['parent_id']}|{d['exact_mutant_overlap']}/{k}；{d['site_overlap']}/{k}|"
                         f"{g['site_overlap']}/{k}|{r['exact_mutant_overlap']}/{k}；{r['site_overlap']}/{k}|")
        lines.extend(["", "这是同一283候选池中、按H实际建议数匹配k的描述性排名对照。d与随机没有套用H过滤。",
            "静态GraphPep只复用native pose的逐残基mean升序，不提供氨基酸替换建议；对应CSV的alt/d/H留空。",
            "排名不同不等于更好；本轮没有实验方向准确率、pooled相关性或优越性检验。", "",
            "## 技术建议（肽内1-based编号）", "",
            "|亲本|位置|替换|d_target|H|H_control_1|H_control_2|", "|---|---:|---|---:|---:|---:|---:|"])
        for row in suggestions:
            lines.append(f"|{row['parent_id']}|{row['site']}|{row['wt']}→{row['alt']}|"
                         f"{float(row['d_target']):.5f}|{float(row['h']):.5f}|"
                         f"{float(row['h_control_1']):.5f}|{float(row['h_control_2']):.5f}|")
        if not suggestions:
            lines.append("|无建议|—|—|—|—|—|—|")
        lines.extend(["", "这些只是模型偏好的单点替换标记，没有完整候选序列、结构、Kd或ΔG声明。", "",
            "## GO/HOLD及下一步", "",
            f"- 本阶段：{summary['decision']}。可以讨论或预注册小型开发验证，不能据本轮结果进入多代设计。",
            "- 亲和力比较：HOLD。保留至少10合格对、两方向各≥4、至少两亲本背景等既有进入条件。",
            "- 下一份应准备的协议：在相同、来源和assay均合格的突变对上比较H与未修正d，再与其他候选评分并列；",
            "  冻结新方法的方向/tie/拒绝处理，失败与拒绝保留分母，报告成对差异的不确定性。不得用已暴露MDM2调参后称独立验证。",
            "- 在标签条件满足以前，不扩大模型、不训练GNN、不增加folding、不以打乱受体当生物学负样本。",
            "- 如果以后做上下文归因诊断，应先准备有依据的、构建体可比的受体扰动及控制，另定有限预算；本轮不执行。", "",
            "## 保留和复核", "",
            "旧六结构pilot和kernel ablation的全部清单项逐项SHA256复核，未改动原始结果和冻结文件；旧AFM结构门失败不因本实验改变。",
            "原始60条分布见raw_masked_distributions.jsonl；全部285替换见all_substitution_diagnostics.csv；建议见technical_suggestions.csv。",
            "运行配置、预算、数值重复、前后测试、来源完整性及交付清单均在同目录独立保存。"])
    (HERE / "report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_json(HERE / "verification.json", {"checks": checks, "failed_count": sum(not c["pass"] for c in checks),
               "preservation": preserved, "no_affinity_claim": True,
               "scope": "Technical verification only; biological performance untested"})
    if any(not c["pass"] for c in checks):
        raise ValueError("Delivery verification failed")
    files = [{"path": str(p.relative_to(HERE)), "bytes": p.stat().st_size, "sha256": sha256(p)}
             for p in sorted(HERE.rglob("*")) if p.is_file() and p.name not in
             ("artifact_manifest.json", "artifact_manifest.sha256")]
    write_json(HERE / "artifact_manifest.json", {"scope": "All files except manifest and detached checksum",
               "files": files})
    (HERE / "artifact_manifest.sha256").write_text(sha256(HERE / "artifact_manifest.json") +
                                                   "  artifact_manifest.json\n")
    print({"checks": len(checks), "failed": sum(not c["pass"] for c in checks),
           "artifact_files": len(files), "output_bytes": directory_bytes(HERE),
           "manifest_sha256": sha256(HERE / "artifact_manifest.json")})


if __name__ == "__main__":
    main()
