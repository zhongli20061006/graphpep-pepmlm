"""Write interpretation only after complete outputs and independent audit pass."""
import csv
import hashlib
import json
from pathlib import Path
import statistics

OUT = Path(__file__).resolve().parent
ANALYSIS = OUT / "analysis"
ROOT = OUT.parents[1]


def read(name):
    with (ANALYSIS / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    audit = json.loads((ANALYSIS / "independent_audit.json").read_text())
    decision = json.loads((ANALYSIS / "decision.json").read_text())
    assert audit["status"] == "passed" and decision["complete"] and decision["technical_pass"]
    original = json.loads((ROOT / "output/foldx_affinity_mdm2_20261006/analysis/decision.json").read_text())
    assert original["continuation_pass"] is False
    rows, structures = read("seed_summary.csv"), read("structure_results.csv")
    lookup = {(r["peptide_id"], r["condition"], int(r["seed"])): r for r in rows}
    assert len(rows) == 8 and len(structures) == 240
    correct = {condition: {pid: sum(lookup[pid, condition, seed]["diagnostic_direction"] == "correct" for seed in (20261011, 20261012))
                           for pid in ("P53_P27A", "P53_L22A")} for condition in ("fixed", "flexible")}
    if correct["flexible"]["P53_P27A"] == 2:
        finding = "释放肽主链后，P27A 在两个种子的三个模型中均恢复了文献预期方向。这支持本方案下 FoldX 对局部结构处理的敏感性。"
    elif correct["flexible"]["P53_P27A"] == 0 and all(lookup["P53_P27A", "flexible", seed]["diagnostic_direction"] == "incorrect" for seed in (20261011, 20261012)):
        finding = "释放肽主链并进行本次局部最小化后，P27A 在两个种子中仍保持与文献相反的评分方向。本方案未纠正这个关键正向对照。"
    else:
        finding = "释放肽主链后，P27A 的方向仍未在两个种子的全部模型中一致恢复；存在混合或种子依赖结果。本方案不足以稳定复现这个关键正向对照。"
    if correct["fixed"]["P53_P27A"] == 2:
        finding += "新固定组也已经方向正确，因此不能将相对旧结果的变化全部归因于主链释放；两臂共同的加氢、侧链处理和力场最小化也发生了变化。"
    if correct["flexible"]["P53_L22A"] < 2:
        finding += "释放组还未在全部种子和模型中保持 L22A 弱结合对照，进一步限制了继续使用该评分目标的依据。"
    forces = [float(r["after_free_force_rms_kJ_mol_nm"]) for r in structures]
    caps = sum(r["iteration_cap_reached"] == "True" for r in structures)
    force_summary = dict(minimum=min(forces), median=statistics.median(forces), maximum=max(forces), threshold=10,
                         accepted=decision["converged_structures"], not_accepted=decision["nonconverged_structures"],
                         reporter_iteration_caps_reached=caps, units="kJ/mol/nm")
    lines = ["# 固定与释放肽主链：结果与解释", "", finding, "",
             f"正式完成 240 个结构、120 个条件内 WT/突变配对和 60 个同一次重复内的 flexible−fixed 差值。48 个预定试跑结果已包含在 240 个中。独立复核通过：{audit['raw_tables']} 张原始表、{audit['raw_numeric_terms']} 个数值字段、所有输入/环境/代码与产物哈希、序列/终端/全精度冻结坐标，以及所有配对聚合。", "",
             "下面是 FoldX AnalyseComplex 的 Interaction Energy 增益，单位 kcal/mol。增益 = 对应 WT reference − mutant；正值表示突变的该代理分数更有利。每个模型先取五个对应配对的中位数，再取三个模型中位数。方向正确要求该种子的三个模型均符合预期，单个汇总中位数不能替代这一要求。", "",
             "|突变|固定：种子 20261011|固定：种子 20261012|释放：种子 20261011|释放：种子 20261012|",
             "|---|---:|---:|---:|---:|"]
    for pid, name in (("P53_P27A", "P27A，预期正增益"), ("P53_L22A", "L22A，预期负增益")):
        values = [lookup[pid, condition, seed] for condition in ("fixed", "flexible") for seed in (20261011, 20261012)]
        lines.append("|" + name + "|" + "|".join(f"{float(r['median_gain_kcal_mol']):+.6f} ({r['diagnostic_direction']})" for r in values) + "|")
    lines += ["", "flexible−fixed 的效果先在同一重复内相减，再进行上述五次和三个模型的中位聚合：", "",
              "|突变|种子 20261011|种子 20261012|", "|---|---:|---:|"]
    for pid, name in (("P53_P27A", "P27A"), ("P53_L22A", "L22A")):
        lines.append("|" + name + "|" + "|".join(f"{float(lookup[pid, 'fixed', seed]['median_flexible_minus_fixed_gain_kcal_mol']):+.6f}" for seed in (20261011, 20261012)) + "|")
    lines += ["", f"力 RMS 验收：{force_summary['accepted']}/240 达到预声明的可动原子向量力 RMS ≤10 kJ/mol/nm；{force_summary['not_accepted']}/240 未达到。结束时该值的范围为 {min(forces):.3f}–{max(forces):.3f}，中位数 {statistics.median(forces):.3f}。报告器记录达到 2,000 次上限的有 {caps} 个。这个单独的较严格验收量与 OpenMM 内部的优化器停止条件有区别，结束或达到迭代上限不能当作通过验收。技术门槛统计结构、坐标和原始评分有效性，验收未达标结果保留并限制科学解释。", "",
              "已有独立 AFM 结构本身已包含 P27A 的末端变化：对齐受体后，B8–12 的 Cα RMSD 为约 2.7–2.9 Å；不能简单归因于 AFM 完全没有预测到构象变化。本次是从同一匹配父代骨架起点出发的局部最小化，允许移动也不保证跨越构象能垒或进入文献关联的螺旋构象。", "",
              "phi/psi 的 alpha-like 区间和 O_i–N_i+4 距离仅为几何描述；Pro 的 N 不作潜在酰胺氢键供体。它们不等于 DSSP、氢键占有率或测得的螺旋比例。五次重复包含重排、加氢与数值处理的变化，不能当作生物学重复或完整热力学系综；图中的范围也不是置信区间。", "",
              "本轮采用 ff14SB/OBC2 隐式溶剂、固定受体主链、弱肽位置约束，以及局部能量最小化。没有游离肽/受体系综、热采样或新 Kd，不能解释为完整 ΔG(binding) 验证；计算一致性也不能证明实验亲和力提升。当前是看过标签后的开发敏感性诊断，原有 continuation_pass=false 保持不变，NCAM1 泛化未验证。", "",
              "结果文件：seed_summary.csv、backbone_summary.csv、pair_results.csv、contrasts.csv；几何记录：torsions.csv、hbond_proxies.csv、residue_displacements.csv；原始数据保存在 ../jobs。计数字段单位为 count，能量字段为 kcal/mol，分项不假定加和等于结合能。", "",
              "文献依据：[Li et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/) 讨论 P27A 的螺旋构象贡献；本轮数据不能单独确认该机制。[OpenMM minimizer API](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html) 说明局部最小化及迭代/容差行为。", ""]
    (ANALYSIS / "interpretation_zh.md").write_text("\n".join(lines), encoding="utf-8")
    summary = dict(finding=finding, direction_correct_seed_counts=correct, force_acceptance=force_summary,
                   original_continuation_pass=False, no_new_Kd=True, no_holdout_validation=True,
                   independent_audit_sha256=hashlib.sha256((ANALYSIS / "independent_audit.json").read_bytes()).hexdigest(),
                   script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ANALYSIS / "interpretation_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
