# Homer1 EVH1 单篇来源筛查（2026-10-20）

结论：**HOLD**。存在同一个 WT 背景的明确双向 BLI 数值，可作为工程六结构 pilot 的资料候选；目前不能按游离短肽启动预测，也不满足完整 benchmark 的至少 4 个增强突变门槛。没有运行模型或评分。

## 来源与附件完整性

原文：Singer, Ramos & Keating, 2024, Protein Science, DOI [10.1002/pro.5094](https://doi.org/10.1002/pro.5094)，[PMC11237882](https://pmc.ncbi.nlm.nih.gov/articles/PMC11237882/)。本次范围仅限这篇文章及其公开附件，未扩展到其他受体或测量来源。

Europe PMC 原始 fullTextXML 已保存。全部五个补充附件经 `supplementaryFiles?includeInlineImage=false` 取得，并与原始 XML 中的五个 MD5 全部核对一致；保存原始 ZIP、XLSX、PDF。S4 只有 EVH1 构建体，S1 为 NGS/引物/PSSM，S2 为 SPOT 设计与信号，S3 为 Dictyostelium motif，补充 PDF 为图 S1–S6。没有找到 BLI SUMO–peptide 完整融合序列或逐重复 Kd。

## 可确认的单点数据

仅使用 Table 2 和 Table 3 中 Homer1 BLI 数值；不混入 ENAH、ddEVH1、coHomer、SPOT intensity。数值为 μM，± 是 **3 次技术重复的均值±SD**；不是 SEM、95% CI 或生物学重复。原文没有说明 WT 与突变是否来自同批独立制备或配对重复。

| 背景及长度 | 原文突变名（肽内一基位置） | WT Kd | 突变 Kd | 仅按点估计 2 倍门初筛 |
|---|---|---:|---:|---|
| DVL3 36 aa | Y(-3)A (Y18A) | 46±6 | 104±8 | 减弱 |
| DVL3 | P(-2)A (P19A) | 46±6 | 145±20 | 减弱 |
| DVL3 | P(-1)A (P20A) | 46±6 | 80±20 | 未过2倍门 |
| DVL3 | P0A (P21A) | 46±6 | 115±10 | 减弱 |
| DVL3 | P1A (P22A) | 46±6 | >150 | 减弱，删失 |
| DVL3 | Y4F (Y25F) | 46±6 | 14±1 | 增强 |
| DVL3 | Y4W (Y25W) | 46±6 | 51±4 | 未过2倍门 |
| DVL3 | Y4A (Y25A) | 46±6 | >150 | 减弱，删失 |
| ABI3 20 aa | G5P (G14P) | 3.9±0.1 | 140±10 | 减弱 |
| DBN1 20 aa | C5P (C14P) | 12±1 | >150 | 减弱，删失 |
| FAT1 20 aa | P5G (P14G) | 98±7 | 15±1 | 增强 |

DVL3 WT 为 `PLPHPGAAPWPMAFPYQYPPPPHPYNPHPGFPELGY`，同背景增强 Y4F 和减弱 P0A 可构成 WT+gain+loss 的数据面候选。全部肽序列和原始 Kd/SD/删失标志保存在 `homer1_kd_records.json`；11 个亲本–单点突变对保存在 `homer1_single_mutants.json`。脚本核验全部为同长度、汉明距离 1。`>150` 保留严格下界，未替换为 150，也未将反向比较另算独立数据。

按原始点估计仅有 **2 gain / 7 loss / 2 uncertain，共 11 个 raw 单点对、4 个背景**。11 不是 eligible 标签数，也不是主指标分母：其中 8 对是有限定量值（2 gain、4 loss、2 未达效应门），另 3 对为 >150 的删失值。即使这 3 个删失值全部合格，剔除 2 个 uncertain 后最多只有 **9 个标签候选**；当前严格已赋标签数为 **0**。这只是来源初筛，并未替父任务生成冻结二分类标签或计算 95% CI；按冻结误差规则审核后数量只可能改变，不能以这些点估计计数宣称过门。完整 benchmark 的至少 10 个 eligible 标签和至少 4 gain 条件均已明确失败；不能把 2 个 uncertain 填入主分母来满足总数门。

## 输入化学和构建体审计

BLI 的肽样本是 **生物素化 SUMO–peptide 融合蛋白**，经 pDW363 在 E. coli 表达、Ni-NTA 纯化，固定在 streptavidin tip。方法 4.3 没有报告切掉 SUMO。没有提供 exact SUMO sequence、融合方向的完整序列、linker、biotin site 或最终融合蛋白末端。因此不能把表中 20/36-mer 的生物学片段静默替换为游离肽。

SPOT 的 N 端 Ac–Gly–Ser、C 端 PEG–cellulose 以及 Cys→Ala 是另一种测量格式，不能套到 BLI。DBN1 BLI 表保留 Cys；S2 SPOT 的 Cys→Ala 不得用于重写 BLI 标签对应的序列。

S4 `EVH1 Domains` 的 Homer1 BLI 表达构建体为 134 aa：

`MHHHHHHSSGVDLGTENLYFQSNAMGEQPIFSTRAHVFQIDPNTKKNWVPTSKHAVTVSYFYDSTRNVYRIISLDGSKAIINSTITPNMTFTKTSQKFGQWADSRANTVYGLGFSSEHHLSKFAEKFQEFKEAA`

方法 4.1 明确 TEV 切割、反向 Ni-NTA 和凝胶过滤。以标准 TEV 位点 ENLYFQ|S 推断切后保留 **SNA + 110 aa Homer1 核心，共 113 aa**；这属于从完整序列作出的合理推断，不是实测 N 端鉴定。S4 核心从 MGEQ 开始、KEAA 结束；本篇范围内未独立核对其 canonical accession 对应范围。没有 MBP 构建体。

## Assay 口径

Octet Red96；室温（没有精确摄氏温度），1000 rpm。buffer 是 EVH1 gel-filtration buffer（20 mM HEPES pH7.6、150 mM NaCl、1 mM DTT）与 PBS buffer（1×PBS pH7.6、1 mM DTT、1% BSA、0.1% Tween-20）1:1 混合。生物素 SUMO-peptide 装载到 0.5–1 nm；每次 association/dissociation 各 100 s，同 tip 依次测 8 个 EVH1 浓度。

用 association 最后10秒均值减 dissociation 最后10秒均值，减去 SUMO-only 对照；Kaleidagraph 非线性最小二乘拟合 y=C[EVH1]/(Kd+[EVH1])。因此是平衡信号拟合 Kd，不是 koff/kon 比。作者特别指出 proline cis/trans 慢相超出该测量时程，信号主要来自游离态中具有适合构象的肽群；不得把这组值解释为已包含充分异构化系综的完整自由能。

## 结构、独立性和进入条件

原文使用的 1DDV 是 Homer1–代表配体结构，不是这些 DVL3/FAT1/ABI3 突变体或 SUMO 融合样本。5ZZ9 是 Homer2–DBN1，受体 paralog 与 BLI Homer1 不一致。此单篇未提供 exact mutant series 的实验结构；未另行宣称完成 PDB 数据库逐链核验。

本项目既往选法暴露、GraphPep 训练集合结构/序列重叠由主任务审计；本筛查不声称盲测或无重叠。原论文做过 motif 或其他模型分析，本身不等于当前项目数据泄漏。

**工程 pilot 数据门有潜力，当前输入门 HOLD**：必须先解决 SUMO/biotin/linker 是否能精确恢复及如何保证模型输入与实验样本对应，再审核标签不确定性与旧开发/训练重叠。完整 benchmark 的至少 10 个 eligible 标签和 4 增强门均失败，不能因其失败就否认双向工程候选，也不能在本篇内拼接 SPOT 或其他 EVH1 受体弥补样本。
