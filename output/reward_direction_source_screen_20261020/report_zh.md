# 候选来源排序与独立方向验证 GO/HOLD 判定

日期：2026-10-20。状态：**HOLD_NO_PILOT**。本轮完成 SKEMPI v2.0 粗筛、ATXN1–Capicua 最高优先级索引命中的原文/结构复核，并与此前完成的 Homer1、PDZ 来源审计合并排序。没有运行新结构预测、模型评分或训练；没有把候选来源称为盲测。

## 当前决定

**暂不启动六结构工程 pilot。** Homer1 是当前最值得保留的外部方向候选：DVL3 同一 WT 背景的 Y25F 增强和 P21A 减弱在已冻结的技术重复误差敏感性规则下通过方向门。但 2024 BLI 实际测量的是固定化的生物素化 SUMO–肽融合物，完整融合序列、连接区和方向未恢复；把游离肽或另一篇论文中的融合前缀作为输入都会改变或假定测量对象。故其输入门不通过。

ATXN1–CIC 是 SKEMPI 粗筛中最像同一短肽/同一 ITC 的候选，但回到原文和行级数据后没有肽增强突变：CIC21 的 W37A、L40S 均为减弱。粗筛中的点估计“增强”来自 ATXN1 受体链单突变；另有双/三突变，不能转作肽突变标签。PDZ 原始补充表有多背景与 15 个肽单突变，但没有一个点估计达到预先冻结的至少 2 倍增强。因此均不能补足双向 pilot 标签。

## 候选来源排序

| 排名 | 来源与证据 | 同受体/assay/WT 方向证据 | 输入与结构 | 判定 |
|---|---|---|---|---|
| 1 | **Homer1 EVH1–DVL3**；Singer et al. 2024，BLI；4 个 WT 背景、11 对原始突变 | DVL3 WT 46±6 μM；Y25F 14±1 μM（r=0.304，增强），P21A 115±10 μM（r=2.50，减弱）。n=3 技术重复；两种冻结的 log-ratio 误差敏感性区间均支持方向。全来源最多 9 个效应大小候选、仅 2 个增强 | 分析物为生物素化 SUMO–肽融合物；2024 完整 BLI 构建体未报告。1DDV 不是该 BLI 样本结构；已知 GraphPep 参考结构/家族暴露，不能称盲测 | **HOLD 输入门**；pilot 暂不 GO。完整面板 NO-GO（数量/类别门不足） |
| 2 | **ATXN1 AXH–Capicua CIC21**；Kim et al. 2013，ITC | 同一 CIC21 WT 下两个肽单点：W37A 8.85e-8→3.82e-6 M（SKEMPI 点值比约 43.2，减弱）；L40S 8.85e-8→1.53e-6 M（约 17.3，减弱）。原文报告 WT 约 100 nM、突变在 μM 范围；未见合格增强。索引中 S602D 的“增强”是 AXH 受体突变 | AXH 表达构建为 ATXN1 aa 563–689、N 端 His tag 经 TEV 切除；ITC 25°C，AXH 25 μM、CIC21 500 μM，100 mM NaCl/20 mM HEPES pH 7.5。CIC21 为 21 aa：EPRSVAVFPWHSLVPFLAPSQ；原文未明确合成肽末端化学。匹配 21-aa 结构为 4J2J；SKEMPI 分组键 4J2L 对应 28-aa CIC28，不能把该结构当作 ITC 肽输入。原文未给重复/误差语义，不能释放严格标签 | **NO-GO**：缺增强；粗筛假平衡由链/突变类型混淆导致 |
| 3 | **a1-syntrophin PDZ–Kv1.5 等肽**；Chen et al. 2008，Supplementary Table 7，FP | 3 个 WT 肽背景、15 对单点。0 个达到 ≥2 倍增强；9 个仅点估计减弱达到阈值、5 个未达、1 个删失。删失不当作负 Kd | 实验 ligand 是 N 端 5(6)-TAMRA–NNG–10 aa core（13 aa）；末端化学/标签接头不能与无标签输入互换。参考 2PDZ 不是这批标记 ligand 的匹配结构 | **NO-GO**：冻结增强阈值下无增强类 |
| 4 | **LGN–FRMPD1 / LGN–FRMPD4**；SKEMPI v2.0 索引命中，SPR | 各 3 行索引记录，点值候选显示 0 增强、2 减弱；尚未逐条回原文确认哪些突变位于短肽链、构建体与误差语义 | 只有链对/索引级信息；未完成输入与结构对应审计 | **低优先级，不进入**：双向门未通过；不为此启动结构任务 |
| 排除 | **MDM2/MDMX–PMI**；SKEMPI 命中 | 同家族/PMI 项目开发面板已暴露，MDM2 还有严格反馈与已知亲和力反例 | 可作开发/压力测试，不是新独立验证集 | **排除独立候选**；不重跑 P27A |
| 排除 | **TCR–pMHC peptide**；SKEMPI 中多组点方向平衡 | 这些条目测量 TCR 对肽-MHC 复合物的结合，不是游离设计肽对目标受体的结合 | 即使 mutation 位于 MHC 呈递肽，也不满足本项目的目标/测量对象 | **排除任务不匹配** |

表中“增强/减弱”始终指 Kd_variant/Kd_WT 的实验方向；GraphPep、FoldX 代理差值不能替代此标签。ATXN1 数字来自 SKEMPI 官方原始 CSV 行级点值；论文原文只明确概述 WT 约 100 nM、两个 CIC21 肽突变在 μM 范围。由于论文没有为此数据提供足以应用冻结误差门的重复/区间语义，这些点值在本项目仍是来源审计证据，不释放为正式方向标签。

## SKEMPI 粗筛的方法和边界

下载并保留官方 `skempi_v2.csv`（SHA-256 见 `verification.json`）。文件 7,085 行，分组后 514 个索引组；按每条有限、精确的 `Affinity_mut (M)/Affinity_wt (M)` 仅做点值 lead 筛选，有 125 组同时出现 ≤0.5 与 ≥2 的点值。此计数**不是**合格的同 assay 短肽 benchmark：SKEMPI 定义为具有复合物结构的蛋白–蛋白相互作用突变数据，组内可混入不同变异链、组合突变和不同分子类型。删失值未补点值；关键词命中仅用于排序，之后必须回查原文和 PDB。

这一轮深入核对最靠前、且表面上最像短肽 ITC 面板的 ATXN1–CIC 后，证实不能通过肽方向门。现有 Homer1 和 PDZ 的详细原始审计仍分别是“输入未解的候选”和“无增强类”，因此当前候选清单没有一个同时达到 pilot 的双向标签、输入/构建体可追溯与测量对应门。

## 沿用的冻结规则与下一步

仍使用 `output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json`，SHA-256：`23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b`。方向按同受体、同 assay、同 WT 背景单点标准氨基酸突变定义；增强要求比值 ≤0.5 且 95% 区间全低于 1，减弱要求比值 ≥2 且区间全高于 1。Pilot 上限仍是 WT+同背景增强+同背景减弱、每者两个固定折叠 seed，共 6 个复合物；seed 沿用 20261008、20261009。主指标仍是体系内、按 WT 背景等权的方向 balanced accuracy；失败、无接触、tie 留在分母，双 seed 均正确才算稳定正确。完整信号门仍为至少 10 对、每类至少 4 对、至少 2 个 WT 背景，并满足冻结的 BA/类别召回/seed 一致率门槛。

**建议的下一步是停止这轮数据库扩展和所有预测计算。** 只有公开资料能独立证明 Homer1 的 2024 BLI 完整融合输入、且能在评分前冻结对应输入表示与命令时，才重新评估其六结构 pilot；否则需要换成构建体、肽末端、同 assay 双向 Kd 与误差都可追溯的新来源。当前不需要训练模型、扩大生成池或重复 MDM2/P27A 计算。PMI/MDM2 的所有结果仍是开发/压力测试，不是盲测。

## 主要来源

- SKEMPI 官方数据库与下载页：[SKEMPI v2.0](https://life.bsc.es/pid/skempi2/)。它描述的是有结构的蛋白–蛋白相互作用突变数据，网页列有 7,085 条 mutation records。
- SKEMPI 2.0 方法论文：[Jankauskaitė et al., Bioinformatics 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6361233/)。
- ATXN1–CIC 原始论文：[Kim et al., Genes & Development 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3613606/)；结构：[4J2J（AXH–CIC21，CIC链21 aa）](https://www.rcsb.org/structure/4J2J)、[4J2L（AXH–CIC28，CIC链28 aa）](https://www.rcsb.org/structure/4J2L)。
- Homer1 来源：[Singer et al., Protein Science 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11237882/)。
- PDZ 来源：[Chen et al., Nature Biotechnology 2008](https://doi.org/10.1038/nbt.1489)。

## 本轮交付

- `candidate_source_ranking.csv`：机器可读排序与状态。
- `decision.json`：当前准入判定和固定计算预算。
- `skempi_index_screen.json`：完整粗筛组级结果，明确是索引 lead 而不是标签。
- `screen_skempi_index.py`：可复核的本地筛查脚本。
- `verification.json`：关键计数、公式重算、原始输入/判据 hash 与本轮零计算核对。
