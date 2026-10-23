# Bcl-w 结构输入与当前预测设置适配性审计

日期：2026-10-20，Asia/Shanghai。状态：**有限资料/输入审计完成；结构失败根因未确定。新预测 0、新放松 0、新 reward 调用 0。**

本目录补充 `PROJECT_HANDOFF_2026-10-20.md`。旧六结构 pilot 已经完成，仍为 **0/6 通过原冻结结构门槛**；不能写成尚未启动，也不能解释为方向准确率为零。原结果、首次失败、冻结配置均保留。

## 判断与本次新增证据

当前可确认：不是已发现的链身份错误，也不是启动命令与 `single_sequence` 产物不符。失败已经出现在原始预测的受体和肽置信度中，不能由后续 GraphPep、FoldX 或 restrained Amber 放松解释。

新增的关键定位是：**受体自身也没有形成可信的预测**。六项受体链平均 pLDDT 为 **39.918–41.324**，六项均没有任何受体残基达到 70；受体内部 PAE 中位数为 **19.69–21.12 Å**。因此当前问题不宜只描述为“肽找不到结合位置”。这些是对已保存 JSON/PDB 的重新解析，不是新结构预测，也不是亲和力评价。

MSA 缺失是较直接的解释候选；加速实现、构建体/天然变体、化学省略与构象异质性尚不能排除。**没有受控新对照，不能确认其中任一项是根因。**

有限决定见 [decision.json](decision.json)：

- **GO：方案准备已完成。** 可以评审一个仅亲本×原两 seed、关闭快速内核路径的新诊断方案；本轮不执行。
- **HOLD：共同受体 MSA 分支。** 六个现有 A3M 没有同源行，不能从中制造一个“共享 MSA”；当前仍缺经审计、可冻结的受体同源 MSA 输入。
- **HOLD：精确实验 pose 分支。** 本次有界结构审查未找到与带标记 assay 分子完全相同的实验复合物。已有相关结构均有明确不匹配。
- **HOLD：完整方向 benchmark、reward 选择和设计闭环恢复。** 标签数量/类别/背景与输入不确定性没有被本次审计解决。

## 1. 实际输入、设置、链映射与置信度

核对范围包括原始 FASTA、六个 A3M、两份 ColabFold 配置与预测日志、执行记录、十二个原始/放松 PDB、六个分数 JSON 及另存的 PAE JSON、放松 manifest 和原 pilot 的 265 文件清单。

| 核查项 | 结果及解释 |
|---|---|
| 原 pilot 完整性 | 265/265 文件尺寸及 SHA-256 匹配；0 偏差 |
| FASTA | 三条受体核心均精确等于本地 UniProt Q92843 1–164；肽为冻结的三条 23-mer；变体分别是亲本位置 13 的 I→L、位置 6 的 I→M |
| A3M | 6/6 为 `#164,23\t1,1`；两条序列化记录分别是受体+23个gap、164个gap+肽；每链各一条生物序列，0 同源序列 |
| 模型/seed/recycle | Multimer v3 model 1；每 seed 三任务；两 seed 为 20261008/20261009；日志均有 recycle 0–3，符合 3 次回收而非提前终止 |
| 生效配置 | 两 seed 共 34 项核对，0 不符；无模板、unpaired、单模型、单seed、fast kernels 与 compile-mode fast 均与冻结设置一致 |
| MSA 上限差异 | CLI 传 256/512，实际配置为 1/1。本机 `batch.py:1363–1368` 在 single-sequence multimer 下自动限缩；属于预期实现行为，不是少执行参数 |
| `host_url` | 默认 API URL 出现在配置中，本身不是调用证据。本机 single-sequence 分支就地生成查询行，现有日志未显示远端 MSA 调用；未做独立网络抓包 |
| 链身份 | 六个原始及六个放松 PDB 均为 A=164 aa受体、B=23 aa肽；序列和连续编号一致 |
| 置信度映射 | 六份 JSON 均为 187 个 pLDDT 与 187×187 PAE；PDB Cα B-factor 与相应 pLDDT 全部一致；另存 PAE 与 score JSON PAE 全部一致 |
| 放松前后 | 原始最小链间重原子距离 0.667–0.985 Å，有严重近接；放松后 2.405–2.606 Å且无<1.5 Å冲突。全复合物 Cα 位移 RMS 为 0.185–0.555 Å；置信度没有重新预测 |

行级证据：[configuration_consistency.csv](configuration_consistency.csv)、[actual_a3m_audit.csv](actual_a3m_audit.csv)、[confidence_and_geometry_audit.csv](confidence_and_geometry_audit.csv)、[recycle_trace.csv](recycle_trace.csv)。本机源码片段及 SHA-256 在 [installed_colabfold_source_excerpt.txt](installed_colabfold_source_excerpt.txt)。

| Seed | 肽 | 受体均值 pLDDT | 肽均值 pLDDT | ipTM（保存JSON） | 原冻结跨链最大 PAE Å |
|---|---|---:|---:|---:|---:|
| 20261008 | parent | 40.685 | 36.443 | 0.12 | 29.78 |
| 20261008 | I3dL | 40.829 | 39.444 | 0.20 | 29.84 |
| 20261008 | I2dM | 40.506 | 36.007 | 0.11 | 29.75 |
| 20261009 | parent | 41.324 | 42.487 | 0.26 | 29.86 |
| 20261009 | I3dL | 40.101 | 41.395 | 0.26 | 30.14 |
| 20261009 | I2dM | 39.918 | 40.817 | 0.20 | 30.28 |

独立重算与旧门槛表六行一致。JSON 的 ipTM 已舍入，日志保留了较多小数，两者不是新的模型结果。原门槛仍为 ipTM≥0.70、肽均值pLDDT≥70、全部双向跨链 PAE 最大值≤10 Å，加身份/接触/放松后冲突要求，**0/6 通过**。

原协议的“界面 PAE 中位数”和本 pilot 的“全部双向跨链最大值”不相同；后者包含远离接触的跨链项，不能偷换为接触局部统计。本次额外给出的全跨链中位数 21.28–23.22 Å、接触局部 PAE 仅作描述，不改变旧门槛。ipTM 与肽 pLDDT 也分别失败，不能靠改 PAE 汇总救结果。

## 2. assay 分子与模型核心代理

完整对应表：[assay_to_input_correspondence.csv](assay_to_input_correspondence.csv)。原文方法支持 N端5/6-fluorescein、C端酰胺的23-mer直接各向异性 Kd；受体区间及 MBP/TEV→c-myc 的来源链由 Dutta 2010 方法核对。当前输入明确没有染料、酰胺、c-myc、接头及其切割残基。来源是 [Foight 2014主文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4168798/)、本地已有补充方法与 [Dutta 2010方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC2896288/)。本次未重下 Dutta 2013 附件。

这里增加一个此前未明确到位点的缺项：**当前模型是 Q133，若干旧结构是 R133；UniProt 把 Q133R 标为 rs910332 天然变体（VAR_048418）**。这不是“旧结构证明模型错”，也不能断言 assay 克隆为 R133。当前 canonical 核心只是可追溯代理；论文方法未提供的完整克隆序列，不能由 PDB 序列代填。[UniProt原始JSON](https://rest.uniprot.org/uniprotkb/Q92843.json) 已留存于 `sources/Q92843_current.json`。

修饰与标签省略会限制向 assay 的解释，但没有带/不带修饰的同条件结构对照，因此不能认定它们造成了全受体低置信度。Dutta 2010 的寡聚体观察属于另一制备，不能移植为 Foight 2014 该批样品的实测比例。FoldX 的 pH/离子强度及标准端基设置也不是原 assay 溶液和化学的完全复现。

来源 ± 仍为平均曲线拟合得到的95%区间；不当作独立 Kd 生物学重复的 SD/SEM。原先的比值敏感性近似不自动释放完整 v2 标签。未标记竞争 Ki 与直接 Kd 继续分开。

## 3. 可用实验结构：相关参考与精确输入分开

RCSB 以 Q92843 accession linkage 的有界查询返回5项；加交接的3FDL，以及原始结构论文指向的4K5A/4K5B，合计审查8个条目。只覆盖该查询与有限引文链，**不是全数据库不存在精确复合物的证明**。新增来源约4.03 MB；既有1ZY3/4CIM/3FDL元数据复用，没有重下模型或软件。

| 结构 | 实际对象 | 本项目可用范围 |
|---|---|---|
| 1ZY3 | Bcl-w相关构建 + BID 20-mer，NMR/分子对接，10个模型；非Bim | BH3沟槽/构象参考；不能当精确assay或纯晶体真值。[数据库](https://www.rcsb.org/structure/1ZY3)、[原论文摘要](https://pubmed.ncbi.nlm.nih.gov/16475813/) |
| 4CIM | 多位点突变Bcl-w + 21 aa突变Bcl-w自身BH3片段，异常二聚/蛋白水解背景，1.50 Å | 不作为本assay的Bim pose。原文支持短链来自自身BH3，且讨论其异常构象。[原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4276867/) |
| 3FDL | Bcl-xL工程构建 + 天然Bim 26-mer，1.78 Å | 仅跨受体的Bim姿态参考；相对A2eT核心为T7A，另有端部长度差异。[数据库](https://www.rcsb.org/structure/3FDL) |
| 1MK3 / 1O0L | 不同突变/尾端的Bcl-w受体NMR结构，无独立肽链 | 受体构象参考，不是Bim复合物。[1MK3](https://www.rcsb.org/structure/1MK3)、[1O0L](https://www.rcsb.org/structure/1O0L) |
| 2Y6W | 含标签/突变的Bcl-w二聚体 | 交换/二聚构象参考，不能把第二受体当肽。[数据库](https://www.rcsb.org/structure/2Y6W) |
| 4K5A / 4K5B | Bcl-w相关构建与DARPins，1.50/1.85 Å | 可结合受体构象参考；非BH3/Bim复合物。原论文称human，当前数据库映射Q1RMX3/Bos taurus；须保留冲突。4K5B entity1 A/B实际是DARPin，Bcl-w在entity2 C/D，不能仅凭名称选链。[4K5A](https://www.rcsb.org/structure/4K5A)、[4K5B](https://www.rcsb.org/structure/4K5B) |

详细结构来源与输入对应表：[structure_source_input_map.csv](structure_source_input_map.csv)。逐残基对应在 [residue_correspondence.csv](residue_correspondence.csv)，完整沉积序列/标签差异在 [deposited_sequence_mapping.csv](deposited_sequence_mapping.csv)。编号统一到**当前输入的1-based位置**：例如1ZY3的沉积注释P116V映射到当前序列P117V，不能把两者当两处突变。

三份已取得 mmCIF 进一步核实：1ZY3第一模型受体/肽 Cα 170/20；3FDL为143/25；4CIM受体 A/B为146/145，短链 P/Q各21。沉积序列长度不等于观察到的坐标长度，缺失残基不能静默补造。4K5A/B本次只审元数据与沉积序列，未下载其坐标。[coordinate_observation_audit.json](coordinate_observation_audit.json)

GraphPep训练列表中的4CIM、Bcl-xL:Bim结构暴露按既有审计披露；不等于本篇Kd标签用于训练。本项目已看标签选体系，仍是开发暴露，不能称盲测。

## 4. 根因假设与可辨别证据

八项假设的支持、反证、缺少的对照与有限检验见 [root_cause_hypotheses_evidence.csv](root_cause_hypotheses_evidence.csv)。优先区分 MSA/结构可识别性与实现路径，再考虑分子表示和构象异质性；不把这些可能性混成一个已确认原因。

历史成功1YCR配置同样是ColabFold1.6.3、fast kernels、BF16/fuse与3 recycle，但使用MSA并有3个参数模型。它削弱“快速内核普遍坏”的说法，也**不能**作为当前Bcl-w同输入的受控阳性。日志中的int64→int32提示和CUDA计时超时提示只说明类型/计时事项，不能凭警告判定结构数值错误。[historical_positive_setup_comparison.json](historical_positive_setup_comparison.json)

没有证据支持低置信度等于实验不结合；没有证据支持重做同样放松能恢复置信度。本次也不重跑240结构主链最小化。

## 5. 已准备、未执行的有限新计算

具体方案：[next_compute_protocol_prepared_v1.json](next_compute_protocol_prepared_v1.json)，配套 [中文说明](next_step_protocol_zh.md) 和 [亲本FASTA](prepared_inputs/parent_only.fasta)。

问题是：固定原核心、两个seed、模型1、3 recycle、single-sequence与无模板，**仅关闭快速内核路径是否足以恢复亲本的可信结构**。使用已保存的同seed亲本作为基线，新增最多2次亲本预测和2个同设置CPU放松结构；不预测突变，不调用任一reward。

预算提案为GPU标记进程累计wall≤30分钟、CPU放松wall≤15分钟、输出≤0.5GiB，串行执行。原结构门槛完全保留，并增加受体均值pLDDT≥70的诊断准入要求。2/2全过只支持这一亲本的实现路径解释；1/2为不稳定；0/2表示关闭该路径单独不足。任何结果都不释放突变评分、完整面板或设计闭环。出现OOM/NaN/程序失败立即停止并保留记录，不换seed/参数救结果。

这项分支优先级来自低成本和可控性，不代表已证实快速内核是主因。MSA分支在缺可冻结输入时保持HOLD；任何“给定实验pose后评分”分支必须另立问题与输入，不能继续称为从序列预测复合物。

本轮授权只覆盖审计与准备，**没有执行该方案**。未来用户明确启动此已准备协议时，按其输入、预算和冻结判据执行，不重复询问已经给出的范围内许可。

## 6. 冻结口径与复核边界

`frozen_direction_criteria_v2.json` SHA-256仍为 `23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b`。正式合格标签的失败、无接触和tie保留分母；当前工程pilot不计算正式方向准确率。GraphPep/FoldX的0.05平局带保留；旧proxy_sign字段并未应用该带，I2dM seed20261009的−0.0093不能被当作稳定减弱。

pilot接触计数为≤6 Å的unique residue pairs，v2 JSON另写distinct heavy-atom pairs；本次两者均列出作为诊断，例如parent seed20261008为53残基对/792重原子对，绝不混成同一baseline。正式reward比较前须另行事前统一，本次不改旧定义。

本轮复核只保证列出的文件、输入对应与数值重算；不保证结构生物学正确，不建立真实ΔG(binding)、新肽Kd、独立benchmark或reward优劣。主报告末尾沿用的Foight2015链接与本轮标签来源Foight2014不同；本目录以2014原文/Table2为准，未改旧报告。

审计规则及授权范围：[audit_scope_and_rules.json](audit_scope_and_rules.json)。来源URL、复用/新取得状态与哈希：[source_retrieval_manifest.json](source_retrieval_manifest.json)。EuropePMC 4CIM论文XML曾返回HTTP500，保留失败记录；随后取得并读取PMC原文HTML，没有把失败响应当证据。

最终文件完整性和交叉检查见 `verification.json`、`artifact_manifest.json`及detached SHA-256。所有写入均限本新目录。
