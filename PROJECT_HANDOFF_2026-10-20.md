# Peptide design 闭环项目：全貌与最新交接

更新时间：2026-10-20（Asia/Shanghai）。本交接由当前工作区的历史交接、来源审计、冻结协议、最新执行记录和原始结果索引交叉核对生成。它补充并更新 PROJECT_HANDOFF_2026-10-07.md，不覆盖旧文件。此次只做状态整理与交接，没有新增研究预测、评分、训练或湿实验。

## 1. 接手时先记住的结论

项目的长期目标是：在约 8 GB 显存、低算力、小规模条件下，验证“PepMLM 生成 → 复合物预测 → 结构/评分与逐残基反馈 → 局部 mask 重生成 → 重新预测评价”的闭环，最终判断反馈是否优于匹配随机，并是否能带来真实实验亲和力改善。不训练大型模型。

当前已经完成一轮反馈重生成的工程链路，但严格 10v10 对照两组代理改善均为 2/10；没有证明反馈引导优于匹配随机。GraphPep 和 FoldX 仍是候选代理，均未通过足以支持闭环恢复的独立亲和力方向验证。新设计肽的实验 Kd、完整 ΔG(binding)、多代优化和 NCAM1 泛化均未验证。

2026-10-20 用户已明确授权并实际完成 Bcl-w/Bim_A2eT 的六结构工程 pilot：6 个原始结构和 6 个放松结构均产生；GraphPep 6 行、FoldX 6 行、ESM2 3 条序列有读数；但 0/6 复合物通过冻结置信度门槛。此轮只能证明流程可执行，不能解释为亲和力方向验证。

当前工作主线是建立可信的 mutation-direction 验证任务，并解决其测量输入与可用复合物结构之间的对应问题。Bcl-w pilot 又暴露了当前结构预测设置的可用性问题，根因尚未定位。不能把项目简化成“换数据集”，也不能宣布必须训练一个新模型。

### 状态速览

| 层次 | 已完成 | 仍未证明/未完成 |
|---|---|---|
| 工程 | PepMLM、AFM/ColabFold、GraphPep、FoldX 和局部 mask 重生成可衔接；一轮重生成完成 | 多代闭环未执行 |
| 反馈策略 | 有探索性 pilot、多父代面板、严格匹配随机对照 | 严格 10v10 未显示引导优势 |
| 评价目标 | 已知 Kd 的 MDM2 开发面板、FoldX 与主链敏感性诊断 | 独立、可比、双向的亲和力方向基准未建成；reward 未验证 |
| 新体系 | 多来源原文/补充表/构建体/重叠审计；Bcl-w 工程 pilot 已执行 | 新 pilot 结构质量 0/6；不能开展有效方向比较 |
| 生物学 | 引用了文献中已有肽的 Kd | 没有本项目新设计肽的亲和力测量；没有湿实验订单 |

## 2. 整个系统中各工具负责什么

1. PepMLM 根据受体和局部 mask 产生候选序列。PPL/语言模型概率表示序列生成指标，不是 Kd。
2. AFM/ColabFold 预测受体–肽相对构象，提供 ipTM、pLDDT、PAE 等结构置信度。置信度高不等于亲和力高，置信度低也不证明实验不结合。
3. GraphPep 对复合物结构评分并输出接触/逐残基反馈。它可以作为结构与接触评价工具；跨序列亲和力目标的有效性尚未建立。接触置信度不是残基结合自由能或突变收益梯度。
4. FoldX AnalyseComplex 的 A/B Interaction Energy 是结合态结构上的能量代理。没有完整游离态体系综、构象熵和热采样，不能称为完整 ΔG(binding) 或精确 ΔΔG。
5. ESM2 masked PLL 只对肽序列逐位 mask 后的原残基 log probability 求和。没有受体条件，只作序列先验候选基线。
6. 界面接触数是几何描述。接触更多不等于更强结合；原子对和残基对计数必须明确区分。
7. 只有在评价目标可信后，反馈选位与匹配随机之间的控制实验才能检验“反馈是否有效”。之后才考虑多代和实验测量。

主线关系：先验证“什么信号可以指导突变”，再恢复“用该信号设计子代”，最后以实验 Kd 检查新设计。结构预测质量与亲和力评价目标是两个分别需要验证的环节。

## 3. 项目如何走到现在

详细历史见 [2026-10-07 交接](E:/Users/graphpep与pepmlm/PROJECT_HANDOFF_2026-10-07.md)。本次重新核对关键报告及生成池快照；没有重跑这些实验。

### 3.1 最初 NCAM1 路线

- PepMLM 产生 1000 条 10-mer，976 条唯一序列；此次再次读取快照核对了行数、唯一数与长度。
- 探索性预测使用 NCAM1 20–718（699 aa）。两个候选各 5 个 AFM 参数模型的平均 ipTM 约 0.370/0.336、肽 pLDDT 20.76/17.76，存在冲突且接口不稳定。
- 同期 1YCR/MDM2–p53 阳性对照可以恢复较可信的复合物结构，证明部分工程链路可工作。
- NCAM1 的结合区域、阳性配体和构建体对应仍不充分，当前暂缓。不能用 MDM2 或 Bcl-w 的工程结果声明 NCAM1 已泛化，也不默认折叠 1000 条全池。

### 3.2 MDM2 闭环与匹配随机

| 试验 | 规模与主要结果 | 当前含义 |
|---|---|---|
| feedback_pilot_mdm2_20261005 | 实验结构反馈引导 3 条 vs 随机 3 条；初始代理改善 2/3 vs 0/3 | 小样本探索 |
| feedback_seed_retest_mdm2_20261005 | 统一放松后旧 seed 3/3 vs 0/3，新 seed 2/3 vs 1/3；有重复序列，采样与折叠 seed 同时改变 | 不算独立新序列复现 |
| feedback_panel_mdm2_20261005 | 3 父代、36 子代，共 39 序列/117 模型；实验结构引导 2/9 vs 随机 2/9，预测结构引导 2/9 vs 随机 1/9 | 未显示可靠普遍优势 |
| feedback_control_1ycr_20261006 | 严格同位点池、突变数、采样预算；引导 10 unique vs 随机 10 unique，60 个子代结构均过门；两组改善均 2/10 | 关键控制实验未证明反馈更优 |

严格对照主汇总是同参数模型配对差值的中位数，不能换成“父代中位数−子代中位数”。原公式的两组平均增益约 +0.01714 vs +0.02084，均只是 GraphPep 代理。

### 3.3 为什么转向亲和力评价目标验证

绕过 PepMLM，用 Li et al. 同一 MDM2(25–109)、同一 SPR 口径的 10 条已知肽，两个家族各 5 条、两 seed×三模型，共 60 个结构。质量门槛均通过后，GraphPep 仍出现局部方向反例：

- pooled GraphPep–log10(Kd) Spearman 约 0.736/0.648；但 p53 家族内部只有 0.100/0.000。
- 七个主要强效突变方向均为 5/7；P27A 实验增强却代理变差，L22A 实验减弱却代理改善。
- 实验坐标参考比较也未解决全部问题，不能把错误一概归因于 AFM，更不能归因于 PepMLM 没有部署好。

随后在同 60 个结构上运行 FoldX：

- seed 20261011：总体 rho 0.976，强效方向 6/7，P27A 错。
- seed 20261012：总体 rho 0.867，强效方向 7/7。
- 联合 continuation_pass=false；P27A 配对代理增益约 −0.114/+0.046 kcal/mol，跨 seed 不稳。总体相关性好不能替代局部突变方向。

L22A 的“成功/失败案例”必须按评分器说明：它是 GraphPep 的方向反例，但在 FoldX 的相应诊断中方向正确。不要笼统归类。

### 3.4 P27A 主链敏感性诊断已经做完

- 固定骨架诊断：24 任务、120 个 WT/突变配对、240 次 AnalyseComplex；P27A 30/30 配对反向。
- 释放肽主链诊断：240 个结构、120 个条件内配对、60 个 flexible−fixed 配对差；48 个试跑已包含在 240 个内。
- 480 张原始表、8160 个数值字段及输入/代码/环境/产物经过独立审计。
- 释放主链仍未稳定恢复 P27A 方向；L22A 对照方向正确。35/240 达到严格力 RMS 验收，205/240 未达。
- 这是 ff14SB/OBC2 下局部最小化，无游离态体系综与完整热采样，不是 MD 自由能验证。当前不重复扩大相同最小化。

## 4. 独立 mutation benchmark 的任务定义

冻结文件：[direction_criteria_v2](E:/Users/graphpep与pepmlm/output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json)，SHA-256：23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b。

第一阶段任务是同一受体构建体、同一 assay、同一亲本背景的肽单点标准氨基酸突变方向分类，先证明增强/减弱方向，不直接做 ΔΔG 回归。

- r = Kd_variant/Kd_parent。增强要求 r≤0.5 且可解释的 95% 不确定区间全部低于 1；减弱要求 r≥2 且区间全部高于 1。
- 误差须说明 SD/SEM/拟合区间、重复数与重复层级。删失只在检测限可证明方向时使用，不填造点 Kd。
- 不跨受体、assay 或亲本拼凑配对；不混 Kd、Ki、IC50、EC50。
- 主指标为体系内、按亲本背景等权的方向 balanced accuracy，另报类别召回、覆盖率、失败与双 seed 稳定性。pooled rho 只能辅报。
- GraphPep 与 FoldX 均按 parent−variant；正值只表示代理改善。v2 的平局带分别为 0.05 score unit 和 0.05 kcal/mol。ESM 按 variant−parent 的 PLL；接触按 variant−parent。
- 正式方向面板的合格实验标签必须完整保留在分母；无接触、低质量、程序失败或 tie 不能被删去以提高准确率。
- 完整单体系信号门：至少 10 个合格配对、两类各至少 4、至少两个亲本背景且类别不由单一背景驱动；BA≥0.70、每类 recall≥0.60、双 seed 同向比例≥0.70。达到也只是有限方法选择信号，之后仍需另选体系确认。
- 六结构 pilot 是一个亲本+一个增强候选+一个减弱候选，各两个 seed 的工程试跑，不是上述性能检验。

原论文用某面板评价其自己的预测器，不自动意味着这些标签进入本项目开发或某个评分器训练。必须分开审计：原论文用途、本项目选法暴露、评分器训练标签/结构/序列重叠、广泛预训练来源。没有 exact 命中不能证明无同源或预训练暴露；看过公开标签选体系后不能称完全盲测。

## 5. 数据来源审计进展

| 来源 | 已核实的价值 | 当前缺口/判定 |
|---|---|---|
| 2021 PMI-L10A/MDM2 | 同 SPR 的 22 个单点 Kd | FP 富集增强候选、开发家族已暴露，误差/精确端基仍有缺口；不作双向独立主基准 |
| a1-syntrophin PDZ | 官方 S7 已取得；18 条肽、3 WT、15 对单点、54 个原始值 | 最优点估计仅约 1.75 倍增强，未达冻结 2 倍门；缺增强类。标签肽实际为 TAMRA–NNG–10-mer，不能用裸10-mer冒充；NO_GO |
| SPSB2–iNOS | 同 assay ITC 与端基信息 | 缺合格增强单点，结构/物种/构建体存在对应缺口；HOLD/不进入 |
| Homer1 EVH1–DVL3 | WT 46±6、Y25F 14±1、P21A 115±10 μM，n=3 技术重复；同背景双向候选 | 冻结的技术误差敏感性区间支持方向，但 BLI 测生物素化 SUMO–肽融合物，2024 完整融合/接头/方向未恢复；输入 HOLD。最多9个效应候选且仅2增强，单篇不足完整面板 |
| ATXN1–CIC21 | 同 ITC 的肽 W37A、L40S 明确减弱 | 索引“增强”来自受体链而非肽；21-mer assay 与28-mer结构分组不能互换；缺增强与误差条件，NO_GO |
| LGN–FRMPD | 索引 lead | 尚无已核实增强与完整原文输入；低优先级，不得当合格数据 |
| Foight BH3/Bcl-w | Bim_A2eT 背景 I3dL、I2dM 有同直接结合 assay 的双向候选；完成有限核心代理输入冻结及6结构运行 | 每受体只有7对等数量/类别/背景不足；拟合CI与化学/标签偏差存在；pilot仅工程用途，现又0/6结构过门 |
| Dutta 2013/Bim_D3fK | 确实存在另一个 Bcl-w 突变背景 | 公开SI没有图4c逐变体直接Kd表与区间，没有合格增强；半衰期表不是Kd表。不能补足完整定量面板 |

PDZ 的“补充表未取得”、Homer1 的“区间方法未冻结”是更早阶段的状态，已被后续工作部分更新；不可再次从零追索已取得附件。Homer1 完整 BLI 融合输入仍缺；Dutta 同一 SI 也已下载读完，不应反复下载充当进展。

GraphPep 重叠审计发现家族/参考结构暴露：Homer1 的 1DDV 在所附 Local62 测试列表，Bcl-xL:Bim 4QVF 在训练集，Bcl-w 4CIM 是突变受体和另一 BH3 肽。它们不等于相应 Kd 标签被用于训练，但必须披露。完整 checkpoint 和 AFM/ESM 预训练来源仍未知。

## 6. 最新 Bcl-w 六结构 pilot：实际完成了什么

主目录：[bh3_bclw_six_complex_pilot_20261020](E:/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020)。[运行结果报告](E:/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020/run_result_report_zh.md)与原始记录优先于同目录的预运行草案/冻结说明。

### 6.1 输入、测量对象和设置

| 角色 | 23 aa 肽序列 | 来源直接 Kd |
|---|---|---|
| 工程亲本 Bim_A2eT | RPEIWITQELRRIGDEFNAYYAR | 39±9.4 nM |
| I3dL 增强候选 | RPEIWITQELRRLGDEFNAYYAR | 12±2.6 nM；r=0.308 |
| I2dM 减弱候选 | RPEIWMTQELRRIGDEFNAYYAR | 260±52 nM；r=6.67 |

来源的 ± 是曲线拟合条件下的95%区间，不是独立 Kd 生物学重复的SD/SEM；比值区间为审计中的事前近似/敏感性估计，不是源论文直接提供的比值CI。正式完整v2标签并未因此自动释放。

assay 用 N-fluorescein/C-amide 肽，以及1–164+c-myc/TEV接头受体；完整接头未恢复。模型使用 Q92843 1–164 核心+标准23-mer，明确省略染料、酰胺及标签，属于核心代理，不是实验分子的精确化学复制。荧光直接Kd不可与未标记竞争Ki混用。

AFM：ColabFold 1.6.3，alphafold2_multimer_v3，model 1，3 recycle，seed 20261008/20261009，每条每seed一个模型，无模板、single_sequence、本地无远端MSA。原始benchmark应用稿曾建议共享受体MSA；此single-sequence更改在新分数出现前记录并冻结，只用于工程pilot。

放松：CPU restrained Amber，max_iterations=200、max_outer_iterations=3、tolerance=2.39、stiffness=10。原始AFM置信度没有因放松而更新。

FoldX：RepairPDB→AnalyseComplex A/B，298 K、pH7.4、ionStrength0.15 M。此标准代理设置与来源assay的pH7.8、50mM NaCl和DMSO等并非完整条件复现，不能把分数单位当实验 ΔG。

### 6.2 结果与判断

- 原始结构6、放松结构6；链序列/编号6/6匹配；接触6/6非零，最小链间重原子距离2.405–2.606 Å，无1.5 Å冲突。
- ipTM 0.11–0.26；肽链平均pLDDT 36.01–42.49；双向界面PAE最大值29.75–30.28 Å。冻结要求分别为≥0.70、≥70、≤10 Å，故0/6通过。
- GraphPep6/6和FoldX6/6产生读数，但保留结构失败标记，仅作为原始执行诊断；没有资格宣称任何配对为稳定正确/改善。
- GraphPep parent−variant：I3dL两seed为−0.0655/−0.1178，I2dM为−0.4672/−0.0093。
- FoldX parent−variant：I3dL为+2.92856/−2.73142 kcal/mol，跨seed变号；I2dM为+5.46553/+1.93856。不能用失败构象支持真实Kd方向或反向机制解释。
- ESM序列PLL的variant−parent：I3dL +0.44469，I2dM −0.39578，仅为序列先验；两个被选候选的符号不能验证reward。
- 受体对齐后跨seed肽CA RMSD：parent3.72 Å、I3dL39.70 Å、I2dM5.37 Å；对应接触Jaccard0.518/0.619/0.392，均为描述指标。

低置信度不证明这些已知文献结合肽“不结合”；也未证明输入有错误。实际失败关口是当前设置下复合物质量，原因尚未定位，不能断言只因single-sequence、标签省略或任何单一模块。MDM2在高质量结构上仍有评分反例，所以解决当前结构问题也不会自动解决reward有效性。

0/6结构通过不等于评分器方向准确率为0，也不允许将六项从将来合格实验标签的正式分母中删除。本次是有限工程pilot，没有计算正式benchmark准确率或pooled相关性。

### 6.3 技术重试、预算与完整性

ESM首次加载因 PyTorch2.7.1默认weights_only=True与fair-esm2.0.0本地loader不兼容而失败，没有产出评分。随后保留原冻结脚本，用独立入口esm_pll_pytorch27.py仅对白名单加入静态检查中唯一的argparse.Namespace；没有weights_only=False。补充记录在第一次成功ESM评分前写入；公式复用原脚本。首次失败和成功重试均在执行记录，测试6/6通过。

实际记录：标记为GPU的任务累计wall time 213.265秒（不是GPU利用率积分或kernel计时），放松与FoldX命令wall time共533.026秒；新增文件约10.75 MB。预算为6结构、2 GPU小时、2 GiB，均未超限。

本次交接再次校验run_result_manifest_v1.json所列265个非临时文件，0个哈希/尺寸偏差。该清单是产物完整性核验；不能等同于旧主链实验的独立数值审计或科研结论验证。

## 7. 当前目标、当前工作与明确下一步

长期目标仍是检验闭环能否优于匹配随机并改善实测结合。当前阶段目标是找到一个测量定义清楚、正反突变充足、结构输入可用的方向验证任务，让候选reward接受检验。当前尚未选出经过独立验证的主reward。

最新pilot完成后，下一个有界任务是“Bcl-w结构输入与当前预测设置适配性审计”，先做只读证据工作：

1. 复核实际FASTA、A3M、ColabFold config/log、链映射及置信度文件是否与冻结命令一致，区分技术错误与模型正常产生低置信度结果。已有链身份检查通过，不能预先宣称映射错。
2. 核对assay受体/标签/肽端基与模型核心表示；明确什么已省略、是否有证据提示会影响结构解释，但不能把未知差异直接当成失败原因。
3. 从原论文和结构数据库审查实验复合物/相关BH3结构：1ZY3、4CIM、3FDL等只能按受体、肽、构建体、突变与方法逐项标为精确/相关参考。不能将BID或另一BH3、Bcl-xL或突变Bcl-w冒充本assay的Bcl-w:Bim结构。
4. 产出结构来源与输入对应表、根因假设/支持与反证表，以及有限下一步GO/HOLD。候选方案可以包括经过审计的共同受体MSA、公开实验结构/模板参考、或转入更适合的体系；它们必须作为新方案冻结，不能事后改本轮门槛救结果。
5. 如果改成在已知实验结合态上评价突变，应明确新问题是“给定pose后的评分”，不能继续称为从序列预测复合物或完整自由能验证。

建议使用新的输出目录，例如 output/structure_input_audit_bclw_20261020；若目录已存在，先核查再新建唯一目录。该审计尚未执行，是下一阶段任务。此次交接未启动新预测。

完整benchmark的数据门也没有消失：Bcl-w没有足够正式双向配对/背景，Homer1输入仍缺，Dutta补件不足。若Bcl-w结构路线不值得继续，应给出证据充分的停止判定，转筛明确缺项的新来源；不无限扩大数据库，不为了凑数混体系。

只有数据/输入/结构可用后，才值得运行预注册的完整方向比较；只有候选reward有有限外部支持后，才回到同位点池、同突变数、同采样预算的引导vs随机设计实验。多代和湿实验在更后面。

## 8. 必须保留的口径差异与风险

- 原v3协议提到界面PAE中位数；本pilot在评分前冻结为所有双向界面项的最大值。不可跨报告混用。当前ipTM和肽pLDDT也各自失败，不能改PAE汇总来宣称pilot通过。
- v2 JSON接触baseline的definition写distinct heavy-atom pairs，协议文字/pilot实际按unique residue pairs计数。未来正式比较前必须统一并事前冻结；保留旧文件，不把两种计数当同一指标。
- pilot的proxy_sign是原始差值正负，并未实施v2的0.05平局带。GraphPep I2dM在seed20261009的−0.0093处于v2平局带，不可直接凭negative字段称稳定减弱。当前所有结构已失败，正式方向分数本来就未计算。
- 原始协议中的GraphPep训练列表“未核验”已被后续限定范围审计部分更新；完整模型预训练独立性仍未知。来源论文用途与本项目泄漏不得混为一谈。
- 没有从柱状图反推精确Kd/CI；半衰期、enrichment、功能活性和Kd分开保存。
- 阈值是工程准入规则，不是适用于所有家族的生物学真值。不能看到结果后调阈值挑成功项。
- 历史各试验按其当时协议解释，不用新的tie/汇总定义回溯改写历史成功率。

## 9. 工作区、环境与运行注意事项

| 用途 | 实际路径 |
|---|---|
| Windows主工作区 | E:\Users\graphpep与pepmlm |
| 同一目录WSL路径 | /mnt/e/Users/graphpep与pepmlm |
| 早期WSL项目/模型目录 | /home/g/pepmlm_pipeline |
| PepMLM权重 | /home/g/pepmlm_pipeline/models/PepMLM-650M |
| GraphPep安装 | /home/g/GraphPep_v1.3 |
| GraphPep/PepMLM Python | /home/g/miniconda3/envs/graphpep_env/bin/python |
| ColabFold Python/CLI | /home/g/miniconda3/envs/colabfold_exp0/bin/python；同env/bin/colabfold_batch |
| FoldX授权二进制 | /mnt/e/Users/graphpep与pepmlm/tools/foldx/foldx5_Linux/foldx_20261231 |
| ESM2 checkpoint | /home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt |
| AFM参数 | /home/g/.cache/colabfold/params/params_model_1_multimer_v3.npz |
| Ubuntu虚拟磁盘 | E:\WSL\Ubuntu\ext4.vhdx（注册表再次核验） |

最近已执行环境：Python3.11.14、torch2.7.1+cu128、fair-esm2.0.0、Biopython1.88；ColabFold1.6.3、JAX/JAXlib0.10.2、OpenMM8.6.1。GPU为RTX5060 Laptop约8151MiB总显存。完整版本和代码/模型SHA-256见pilot冻结配置与launch_preflight_v1.json。

注意 /mnt/e 工作区与 /home/g/pepmlm_pipeline 不是同一副本。旧副本feedback_panel.py哈希不同，且缺部分新脚本；没有同步。当前pilot明确调用/mnt/e下pipeline，并指定输出目录内resmap兼容文件。不要为了读报告运行scripts/sync_pipeline.sh或盲目同步。

本次交接时Windows快照：C盘约4.17GiB空闲，E盘约109.95GiB空闲。此前C盘只剩277MiB的说法是更早快照，已发生变化。WSL的df显示虚拟文件系统容量不等于宿主盘实际剩余；本pilot输出与Ubuntu虚拟盘均在E。未来计算启动时重新检查，不沿用任一快照当资源保证。

E:\AGENTS.md是FastAPI专项规范，此研究目录没有对应FastAPI应用；本次未发现工作区内AGENTS.md。新会话仍需检查当时的适用上级/子目录规则。

## 10. 证据阅读顺序与旧状态处理

1. 本文件，再读最新pilot结果与机器记录：
   - [run_result_report_zh.md](E:/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020/run_result_report_zh.md)
   - [structure_gate_audit_v1.csv](E:/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020/structure_gate_audit_v1.csv)
   - [paired_proxy_summary_v1.json](E:/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020/paired_proxy_summary_v1.json)
   - [run_execution_record_v1.json](E:/Users/graphpep与pepmlm/output/bh3_bclw_six_complex_pilot_20261020/run_execution_record_v1.json)
   - 同目录run_config_frozen_v1.json、launch_preflight_v1.json、run_result_manifest_v1.json及sha256、ESM补充记录、原始seed/score/job输出。
2. 来源与输入：
   - [Foight/Bcl-w来源审计](E:/Users/graphpep与pepmlm/output/bh3_bclw_source_audit_20261020/report_zh.md)与input_representation_frozen_pilot_v1.json。
   - [Dutta补充材料最终追查](E:/Users/graphpep与pepmlm/output/bh3_second_background_audit_20261020/supplement_followup_zh.md)。
   - [候选来源排序](E:/Users/graphpep与pepmlm/output/reward_direction_source_screen_20261020/report_zh.md)。
   - [Homer1/PDZ审计](E:/Users/graphpep与pepmlm/output/independent_direction_audit_20261020/report_zh.md)与followup_20261020/report_zh.md、independence_addendum_zh.md。
3. 冻结方向方案：output/reward_direction_benchmark_20261008/protocol_draft_zh.md与frozen_direction_criteria_v2.json。
4. 解释历史瓶颈的四份报告：
   - [严格反馈控制](E:/Users/graphpep与pepmlm/output/feedback_control_1ycr_20261006/report_zh.md)
   - [GraphPep亲和力诊断](E:/Users/graphpep与pepmlm/output/affinity_diagnostic_mdm2_20261006/report_zh.md)
   - [FoldX亲和力开发诊断](E:/Users/graphpep与pepmlm/output/foldx_affinity_mdm2_20261006/analysis/report_zh.md)
   - [主链敏感性解释](E:/Users/graphpep与pepmlm/output/backbone_diagnostic_20261007/analysis/interpretation_zh.md)

旧HOLD_NO_PILOT、尚未评分、尚未下载等字段按生成时点解释。用户后来已授权并完成这一Bcl-w工程pilot；不能因旧decision把当前说成零预测，也不能把一次有限授权扩大成新试验/全池授权。run_config_frozen_v1.json中的NOT_STARTED是原始冻结状态，保留不改；实际执行以run_execution_record_v1.json和结果为准。

旧stage1_report、Experiment0 README、FoldX README存在已知过时状态，按最新结果与原始证据判断。最新人工报告也不能覆盖冻结协议、原始CSV/JSON和完整性限制。

## 11. 授权、偏好与交接边界

- 中文沟通，准确区分已完成、建议、已验证、仍未知。低成本、可逆的资料/输入核查和文档整理直接推进，不反复要求确认。
- 用户此前已授权并完成本次6结构pilot；不重问它的历史许可，也不重跑已完成单元。新计算方案应先做成可评审的具体协议与输入；用户若明确启动，按既有授权执行，无需重复确认。
- 此次请求是项目全貌与新会话交接，没有授权另起一个结构试验。新会话提示词授权继续有限资料/输入审计，默认不跑新预测。
- 不重装/重下GraphPep、FoldX或已有模型，不默认云算力、付费、第三方联系、湿实验下单。
- 新实验用新目录，保留原始结果/首次失败/冻结配置。不要修改旧结果来使manifest重新通过；不要重跑全pipeline仅为重新生成报告。
- 不假设旧代理、下载或后台任务仍存活；需要时只读核对实际进程。
- 记忆中的旧状态只能作线索，本文件与当前机器记录优先。没有获准修改长期memory，本次仅写工作区交接文件。

## 12. 本次交接产物与核验范围

新建PROJECT_HANDOFF_2026-10-20.md和NEW_SESSION_PROMPT_2026-10-20.txt；旧2026-10-07文件保留。本次核验包括：关键文件存在、最新265文件manifest哈希、NCAM1快照1000/976计数、pilot机器记录/报告与数值、C/E容量和Ubuntu存储位置。没有重新生成分子、重新折叠/评分，也没有为所有历史试验重新执行独立数值审计。

交接文件与引用来源的完整性快照见 output/project_handoff_20261020/verification.json。新会话应先复核本文件，再完成第7节的有限结构输入审计并形成GO/HOLD，而不是从零重新筛一遍已排除文献。
