# PMI 补证审计：开发面板用途与边界

日期：2026-10-08  
状态：只做原文与现有输入审计；没有为该面板运行新结构预测或评分。  

## 审计结论

2021 年 Li 等人的 PMI-L10A 面板可以作为既有 MDM2/PMI 家族上的**增益侧流程检查**，不构成新受体体系或独立双向方向基准。论文先用 FP 从 94 个 PMI-L10A 单点类似物中筛得 22 个比 PMI-L10A 结合更强的类似物，再用 SPR 报告这些 22 条及母肽的 Kd。因此 Table 2 是按增强富集的 Kd 面板，未测得的其余类似物不能反向充当减弱标签。该论文同时明确将原始 PMI 的 Phe3、Tyr6、Trp7、Leu10 作为关键位点，并对 Leu10 做 Ala 替换以降低亲和力，方便后续定量；这仍属于同一 MDM2-PMI 开发家族。

Table 2 中 MDM2 的 PMI-L10A Kd 为 897 ± 33 nM；代表的单点变体包括 T1L 123 ± 10 nM、A4L 98 ± 5 nM、N8A 229 ± 11 nM、S11M 198 ± 4 nM 和 P12Q 132 ± 2 nM。论文的 SPR 方法为 25°C、10 mM HEPES、150 mM NaCl、0.005% P20、pH 7.4，使用 MDM2(25–109)；每项三个独立实验、每个实验重复两次。可访问正文没有定义表格 ± 是 SD、SEM 还是拟合误差，因此不能把这些 ± 直接当作冻结规则所需的 95% 配对不确定区间。

该文肽以 Boc 固相法合成并经 HF 切割，但未说明每条 SPR 分析物精确的 N/C 端封闭状态或树脂化学。用于竞争 SPR 的检测肽是溶液中的未标记竞争物；TAMRA 标记用于 FP tracer，不能误记到 SPR 分析物上。5UMM 是 MDM2–PMI-M3 复合物，不是 PMI-L10A；3EQS 与 3LNZ 分别提供母本 PMI 和 N8A 的相关结构，但不等同于 L10A 及 Table 2 变体的实测复合物。

## 与当前项目的关系

工作区此前已有原始 PMI/MDM2 和 p53/MDM2 的 Kd 面板，并据其检查/选择 GraphPep 亲和力评分与反馈规则。2021 年 PMI-L10A 的具体序列在现有结果表中未见评分记录，但受体、肽家族与项目选法相同，故必须标为 `DEVELOPMENT / WITHIN-FAMILY STRESS TEST`。不能通过将父代与增强子代反向排列，把一个方向富集的发现集伪装成独立的双向样本。

现阶段允许用途：在未来确实补齐端基与误差定义后，对冻结工作流做公开开发面板上的增益方向敏感性检查，并与历史 GraphPep/FoldX 结果并列展示。禁止用途：新体系独立验证、盲测、证明普遍亲和力能力或据此选择主评分器。

## 进入门仍未满足

- 双向证据：22 条 SPR Kd 变体均属于经 FP 选出的增强集合；没有可确认的同 assay、同父代减弱 Kd 面板。
- 不确定度：论文报告独立实验与技术重复数，但没有说明 Table 2 的 ± 统计语义或配对结构。
- 结构/输入：Table 2 变体无同序列实测复合物结构，端基状态未从合成方法明确恢复。
- 独立性：同 MDM2/PMI 家族已参与当前项目方法选择；GraphPep WSL 训练集重叠审计仍不可读，不能作无重叠声明。

因此，PMI 补证到此停在来源层面，不解冻预测或评分。下一项高价值工作是找出一个不同受体家族、同一 Kd assay 中真正同时有增强与减弱单点变体、端基和构建体可追溯的面板；如果找不到，就应正式记录当前公开资料不足，而不是扩大计算。

## 原始来源

- Li et al. 2021，方法、Table 2/3、PDB 5UMM/5UML：[原文 PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463443/)
- 原始 PMI/Mdm2-PMI-N8A 开发面板：[Li et al. 2010 原文 PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/)
- 相关结构：[3EQS](https://www.rcsb.org/structure/3EQS)、[3LNZ](https://www.rcsb.org/structure/3LNZ)、[5UMM](https://www.rcsb.org/structure/5UMM)
- 本地已存在的逐肽输入：`output/independent_affinity_validation_20261007/candidate_affinity.csv`、`mutation_pairs.csv`、`candidate_source_audit.csv`；旧目录未修改。
