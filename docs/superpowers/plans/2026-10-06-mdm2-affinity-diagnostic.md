# MDM2 实测亲和力与结构来源诊断 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 绕过 PepMLM，依次检验 AFM→GraphPep 的跨序列亲和力排序、实验结构对照，以及逐残基选位规则是否值得继续。

**Architecture:** 使用独立输出目录，复用现有 AFM、Amber 和 GraphPep 工具链。实测数据只取 Li et al. 2010 同一 SPR 竞争实验、同一 synthetic MDM2(25–109) 构建体的十条 12-aa L 肽；先冻结数据与判据，再运行预测和评分。实验结构比较仅使用化学序列匹配的 PMI/PMI-N8A，并显式处理缺失末端和晶体多副本。

**Tech Stack:** 现有 WSL graphpep_env、colabfold_exp0、BioPython、pandas、SciPy、matplotlib；本地 RTX 5060 8GB；不训练模型。

---

## Task 1: 冻结实验基准

**Files:**
- Create: `output/affinity_diagnostic_mdm2_20261006/reference_affinity.csv`
- Create: `output/affinity_diagnostic_mdm2_20261006/protocol.json`
- Create: `pipeline/affinity_diagnostic.py`

- [x] 独立核对 https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/ 的 Table 1/3：PMI-WT、N8A、Y6A、F3A、W7A；p53(17–28)-WT、P27A、K24A、L22A、L26A。原始 Kd/SD 从 M 转 nM，保存原表与 assay 信息；不混用 ITC、MDMX 或不同长度数据。
- [x] 记录 PMI `TSFAEYWNLLSP` 与 p53 `ETFSDLWKLLPE` 的真实单点序列；使用共同 85-aa MDM2 序列。读取所有原始 PDB 结构副本，不根据 GraphPep 分数挑选。
- [x] 在评分前冻结 reference CSV、脚本、MSA、GraphPep 源码哈希与判据。主终点为 GraphPep 对 `log10(Kd_nM)` 的正向 Spearman；GraphPep 越低、Kd 越低表示方向一致。

## Task 2: 绕过生成器的 AFM→GraphPep 基准

**Files:**
- Create: `scripts/run_affinity_diagnostic.sh`
- Create: `output/affinity_diagnostic_mdm2_20261006/afm_msa/`
- Create: `output/affinity_diagnostic_mdm2_20261006/analysis/`
- Test: `tests/test_affinity_diagnostic.py`

- [x] 使用共同受体 MSA、AFM multimer-v3 模型 1/2/3、seed 20261011 与 20261012、3 recycles、不使用模板，预测十条已知肽，共 60 个模型。所有模型统一 Amber 有约束最小化。
- [x] 质量门槛与旧轮相同：ipTM≥0.7、肽 pLDDT≥70、界面 PAE≤10 Å、最小界面重原子距离≥1.5 Å、有接触。保留未通过模型的评分用于诊断，主质量过滤结果要求每个 seed 至少 2/3 可用。
- [x] 对总体、PMI 家族和 p53 家族分别计算排名；比较两个 seed 的稳定性。不得把六个模型当作六条独立肽或生物学重复。
- [x] 冻结七项有≥3倍实测变化的方向检查：N8A-PMI、P27A-p53 应变强；Y6A/F3A/W7A-PMI 与 L22A/L26A-p53 应变弱。K24A 小于 3 倍，仅描述性报告。
- [x] 继续亲和力闭环的操作性门槛：每个 seed 的总体 rho≥0.6、两个家族 rho 均>0；七项方向检查至少六项正确，且两项变强对照均正确、所有检查可比较。此门槛不对应统计显著性或普遍适用性。结果：两个 seed 均未通过。
- [x] 用实际数据独立核验 nM 单位、评分方向、逐模型配对、模型数、原始接触公式、冻结哈希。添加只覆盖上述科学计算错误风险的少量测试。

## Task 3: 实验结构与 AFM 结构对照

**Files:**
- Create: `output/affinity_diagnostic_mdm2_20261006/pdb/3LNZ.cif`
- Create: `output/affinity_diagnostic_mdm2_20261006/native_complexes/`
- Create: `output/affinity_diagnostic_mdm2_20261006/analysis/native_vs_predicted.csv`

- [x] 复用 3EQS；下载 RCSB 3LNZ。按几何距离把各肽链配到受体链，验证化学序列与构建体，保留所有副本。
- [x] 原始实验重原子坐标不移动，标准化链 A/B 和残基编号；记录缺失肽末端和受体位置，不人为补坐标。
- [x] 主结构来源对照固定为所有晶体副本共同已解析 MDM2 27–108（82 aa）和 PMI 肽位点 1–11；3LNZ 部分副本缺失 25/26/109，不能用 26–109 范围保留全部副本。AFM 预测来源于完整 85+12-aa 输入，再按相同位置裁剪后评分。PDB 原样解析范围、完整 AFM 和裁剪 AFM 结果均保留。裁剪评分仅为共同坐标比较，不是对应截短肽的实测 Kd。
- [x] 比较 PMI 与 N8A 的方向和晶体副本分布，并计算受体对齐后的肽 CA RMSD。仅两个不同序列，不能外推普遍排序准确率。
- [x] 1YCR 化学肽与 p53(17–28) 不同，不能给其 native score 赋予本基准的 440 nM；仅作为已有构象参考，排除精准亲和力归因。

## Task 4: 选位规则诊断与交付

**Files:**
- Create: `output/affinity_diagnostic_mdm2_20261006/report_zh.md`
- Create: `output/affinity_diagnostic_mdm2_20261006/analysis/feedback_vs_mutation.csv`
- Create: `output/affinity_diagnostic_mdm2_20261006/analysis/affinity_scatter.png`

- [x] 在已知 parent 的质量通过 AFM 模型上汇总逐残基接触置信度、接触数与缺失反馈位置；与实测单点突变方向对照，检验“低接触置信度优先突变”是否漏掉 N8A 等非接触位置。
- [x] 若基准没有通过预定门槛，停止将当前 GraphPep 指标作为亲和力迭代目标，完成诊断报告与后续独立目标方案；不通过换 seed、改指标或重新选样宣称成功。
- [x] 已检查进入下一轮的条件：上述门槛未通过，因此未执行新子代生成或选位/替换策略优化；本轮不以小样本基准宣称 ΔG 或新子代实测亲和力提升。

## Run commands

```bash
/home/g/miniconda3/envs/graphpep_env/bin/python tests/test_affinity_diagnostic.py
bash scripts/run_affinity_diagnostic.sh
```

运行产物须有十条唯一已知肽、六十个 AFM 模型、每个模型的质量与 GraphPep 状态、冻结协议与核验文件、原始实验结构副本及对照表。仓库无 Git 元数据，使用新输出目录隔离，不创建提交。
