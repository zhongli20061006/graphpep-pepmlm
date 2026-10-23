# Bcl-w 第二突变背景来源补审

日期：2026-10-20。范围为只读文献/输入审计；未运行预测、评分、训练或实验。

## 判定

Dutta 等 2013（PMCID PMC3631442）确实提供了 Bcl-w 对 Bim BH3 单突变系列的直接荧光偏振/各向异性 Kd 测量，构成 Foight 2014 之外的一个工程化亲本背景 `Bim_D3fK`。它补上了“第二背景存在”的证据，但**没有补齐可进入冻结方向 benchmark 的数据条件**：论文图 4c 中该背景的三个变体只有图形化 Kd，无可提取的逐变体数值/区间/误差；正文只明确给出 E3gV 对 Bcl-w 约 2.5 倍减弱。F4aL 的文字支持对 Bcl-w 的相对减弱但没有定量幅度，I3dA 被描述为多个受体（含 Bcl-w）可容忍，不能强行判方向。因此此来源不能补足增强类，也不能达到每背景定量、可复核标签的标准。

**当前结论：完整独立方向 benchmark 仍为 HOLD；这不是项目停摆。** 可保留 Foight Bcl-w/Bim_A2eT 的小型工程可行性 pilot 配置供后续决定，但本轮用户授权仅覆盖来源审计，故没有运行。Dutta 数据不能作为新一轮盲测：它的研究设计明确使用既有 Bim BH3 SPOT 突变数据/PSSM 指导筛选和选取突变；本文新测的直接 Kd 虽不等于旧 SPOT 标签，但存在明显开发/选择暴露。

## 数据来源审计

| 项目 | 审计结果 | 对 benchmark 的影响 |
|---|---|---|
| 论文/数据 | Dutta et al., ACS Chem Biol 2013, DOI 10.1021/cb300679a, PMC3631442；图4c测 Bim-D3fK 背景的 I3dA、E3gV、F4aL | 同一 Bcl-w受体上的突变子系列；仅作补充来源 |
| assay | 荧光肽与纯化受体直接结合，single-site拟合；论文称 direct-binding；37°C，20 mM phosphate、50 mM NaCl、1 mM EDTA、0.001% Triton X、5% DMSO、pH 7.8 | assay 类型相容，但条件不同于 Foight；不可简单合并绝对Kd |
| 肽长度/端基 | Bim BH3 variants 为23 aa；直接检测用荧光标记肽。该文明确FP competition用 N-acetyl/C-amide，但图4c直接标记肽的末端处理未在可见方法句中逐项明确；需以补充材料/原始合成表确认 | 端基/染料细节未完全核实，不能作为严格匹配标签 |
| 受体构建 | 采用该组前序研究相同的pro-survival protein constructs；本篇未在本段重列Bcl-w精确残基和标签；引用Foight等方法 | Bcl-w精确构建体与Foight的等同性未由本篇独立证明 |
| Kd与误差 | Figure 4c仅绘制各突变Kd；正文仅明确 E3gV 对 Bcl-w 约2.5倍减弱。未见图4c变体绝对Kd、置信区间或误差定义；不可从柱高反推精确数值/不确定度 | 可记录一个定性/文本方向，不纳入冻结的定量主分析 |
| 双向证据 | E3gV为减弱；F4aL正文描述其相对Bcl-w的特异性降低，但无明确折叠数；I3dA在Bcl-w等受体被描述为可容忍。没有合格的增强标签 | 不满足增强+减弱平衡 |
| PDB | 文章结构讨论引用Bcl-w结构1O0L，并将其与Bfl-1:Bim结构2VM6对齐；没有报告Bcl-w:Bim-D3fK复合物晶体结构 | 不能把结构参考误称为实测突变复合物 |
| 开发/选择暴露 | Bfl-1设计流程明确从已发表的Bim BH3 SPOT突变数据/PSSM出发；本文随后对工程筛选序列做突变拆分。公开数据已被用于选择/解释突变 | 非盲测，且存在选择偏倚；不做模型外部验证的独立证明 |

图4c原图已保存在 `raw/Dutta2013_Figure4.jpg`，SHA-256：`902fe9488e84e52dc150f84d9ed6da101523fef912a478c7f381f60a9d2a9fec`。审计只从正文抽取明确报告值；未对图形柱高数字化。

## 更新后的决策和后续门槛

1. 保持已冻结的方向判定：增强 `Kd_variant/Kd_parent ≤ 0.5` 且比值95%区间完全低于1；减弱 `≥2` 且区间完全高于1；其余或误差/构建体/末端不清时标 `HOLD`。
2. 完整 benchmark 仍需至少10个合格配对、增强/减弱各至少4、至少2个亲本背景，且两类不能由单一背景驱动。Foight目前每个受体仅7对；Dutta图4c现有公开抽取不能增加符合定量门槛的配对。
3. 资料补齐的具体缺项：Dutta图4c各变体的逐重复曲线/数值、拟合Kd及误差/CI；图4c荧光标记和肽末端完整化学；精确Bcl-w构建序列与是否与Foight一致。若supplement没有这些材料，停止把Dutta纳入定量主分析。
4. 因此暂不生成新的六结构集合，也不跑预测/评分。先完成针对Foight已冻结 pilot 配置的命令、软件/模型权重哈希、GPU/时间预算与无接触失败规则冻结；完成后仍须把“工程输入可重复性 pilot”和“affinity方向验证”分开报告。

## 可引用原始来源

- Dutta et al. 2013 全文及方法/数据：<https://pmc.ncbi.nlm.nih.gov/articles/PMC3631442/>；Figure 4c原图：<https://cdn.ncbi.nlm.nih.gov/pmc/blobs/e9f7/3631442/7aff725b1081/nihms-443232-f0004.jpg>
- Foight et al. 2015 Bcl-w/Bim_A2eT直接Kd对照和材料方法：<https://pmc.ncbi.nlm.nih.gov/articles/PMC4520770/>
