# FoldX MDM2 亲和力指标验证

**正式 FoldX 5.1 评分已完成：60/60 个模型成功，0 个失败。** 固定的 affinity screen 门槛未通过：seed 20261012 通过，seed 20261011 未通过；独立原始输出与统计复算通过。结果见 [完整报告](analysis/report_zh.md)、[决策 JSON](analysis/decision.json) 和 [独立评分审计](analysis/independent_score_verification.json)。

60 个模型的源结构哈希、导出重原子坐标和完整 85+12-aa 序列均通过 [独立输入核验](independent_input_verification.json)。FoldX 原始表格均成功解析，并由独立脚本从原始输出重算 40 条汇总、32 条突变方向、12 组排名相关和两个 seed 门槛。8 项技术与科学计算测试通过；测试中的合成能量只存在于临时目录，不属于真实结果。

`runtime_preflight/` 是正式批处理前用一个单独输入副本完成的 FoldX 5.1 运行检查。那里生成的文件和试算不属于 60 个正式模型，也不参与任何相关系数、配对增益或通过判定。

输入复用十条有同一 assay 实测 Kd 的 12-mer 肽及已完成的 60 个 AFM 复合物模型：两个 seed × 三个模型。原始数据位于 `../affinity_diagnostic_mdm2_20261006/`。

## 固定计算流程

1. 核对原始结构哈希，导出不移动坐标的重原子输入；受体 A 链 85 aa，肽 B 链 12 aa。
2. 每个模型独立 RepairPDB；验证修复后化学序列仍完整，再 AnalyseComplex(A,B)。
3. 统一 temperature=298 K、pH=7.4、ionStrength=0.15 M、water=-IGNORE、noNterm=false、noCterm=false，保持游离末端处理。298 K 是对文献室温的统一近似；0.15 M 代表主要 NaCl 浓度，未精确复现全部 assay 缓冲液。
4. 主指标为 Interaction Energy，单位 kcal/mol；配对增益 = 同 seed、同模型号的父代能量 − 突变能量，取差值中位数。

不将 BuildModel 的稳定性差值当作结合能；不将 FoldX 能量值直接宣称为精确实验 ΔG/ΔΔG。这里先检查方向和排名是否能支持优化。

## 评分前固定的门槛

每个 seed 十条肽全部可比较、总体 Spearman(energy, log10Kd)≥0.6、PMI/p53 各自>0；七项有至少三倍实测变化的突变中至少 6/7 方向正确，N8A 与 P27A 两项增强必须正确，全部七项须可比较。每个 seed 的质量过滤统计及配对至少需两个有效模型。

实际结果：seed 20261011 的 P27A 方向错误（配对 gain = −0.114 kcal/mol），该 seed 为 6/7 且未满足 P27A 必须正确的门槛；seed 20261012 的七项方向全部正确，但 P27A gain 仅 +0.046 kcal/mol。总体继续判据为 FAIL；这不构成精确 ΔG 预测，也不据此宣称 holdout 验证。

固定 60 个模型必须全部完成评分尝试后才作继续决策；运行中快照即使暂时满足排名或方向门槛，也不能提前宣布通过。

|关键对照|父代 Kd (nM)|突变 Kd (nM)|预期方向|
|---|---:|---:|---|
|PMI N8A|3.2|0.49|增强|
|p53 P27A|440|51|增强|
|p53 L22A|440|5000|减弱|

来源：[Li et al. 2010 Tables 1/3](https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/)。既有十条数据是方法开发/筛查面板，不是盲测集，也尚未审计模型训练集重叠。

## 条件性后续

只有新目标通过筛查，才进入未参与选方法或调参的额外已知突变验证；通过后再恢复小规模引导/随机子代比较，最终以同一条件下实测 Kd 检验。此处尚未执行新子代实验或 MD/MMGBSA，不宣称新子代亲和力或 ΔG 改善。

运行方式见 [tools/foldx/README.md](../../tools/foldx/README.md)，协议与进度见 [实施计划](../../docs/superpowers/plans/2026-10-06-foldx-affinity-validation.md)。
