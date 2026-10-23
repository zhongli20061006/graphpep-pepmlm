# 随Git纳入的阶段数据

日期：2026-10-23。用户明确要求代码之外也纳入现阶段的一些本地数据。

本次使用目录白名单，保留支持、失败、拒绝、HOLD和负面结论；不按结果优劣挑选行。原工作区其他文件不删除、不移动。逐文件清单、尺寸、SHA256及目录统计见同目录`repository_data_manifest.json`，所列路径均相对仓库根目录。

初始快照共808个数据文件，34,715,464字节，约33.1 MiB；代码和说明文档另外计入Git。

## 1. 近期主要执行记录

| 目录（位于output） | 纳入内容 | 结论边界 |
|---|---|---|
| `algorithm_hypothesis_plan_20261022` | H初始提案、来源/输入、控制和资产哈希记录 | 提案中的NOT_RUN按生成时点解释 |
| `receptor_context_mutation_probe_v0_20261022` | 原始masked概率分布JSONL、建议/基线表、冻结配置、数值核对、日志、测试和代码 | 技术原型；不是新设计肽的亲和力证据 |
| `receptor_context_pmi_development_comparison_20261022` | 新/复用分布、七对原始比较、来源审计、误差敏感性、预算和失败记录 | 七对点效应开发比较；正式标签未释放 |
| `receptor_context_disagreement_audit_20261022` | 361项总表、84项分歧、精确/近似实验匹配、谱系和重复来源判断 | 后验开发审计；没有合格分歧标签 |
| `bh3_bclw_six_complex_pilot_20261020` | 六结构原始/放松PDB、置信度JSON、A3M、GraphPep/FoldX/ESM2读数、日志和冻结文件 | 0/6结构过门；不计算正式方向准确率 |
| `bclw_parent_kernel_ablation_v1` | 亲本双seed结构/置信度、放松结果、预算对账和完整性记录 | 已完成，0/2过门 |
| `bclw_construct_msa_audit_20261022` | 原检索元数据、93条UniProt候选TSV、事前筛选规则和来源摘录 | 检索当时未完成比对，已由次日审计续接 |
| `bclw_msa_input_audit_20261023` | 逐候选去向、全长比对证据、核心映射、22条代表、A3M、多样性、比对歧义、输入门和未来协议 | 输入GO_PREPARE_ONLY；预测未启动 |
| `structure_input_audit_bclw_20261020` | 原FASTA/A3M/配置/链映射审计、UniProt/RCSB元数据及1ZY3/3FDL/4CIM坐标 | 相关结构不是精确assay结构；完整论文HTML未纳入 |
| `reward_direction_benchmark_20261008` | 冻结方向标准、协议草案及来源审计 | 不事后改门槛 |
| `reward_direction_source_screen_20261020` | SKEMPI官方CSV、候选排序及来源筛查记录 | 数据库索引命中不等于合格短肽标签 |
| `project_handoff_20261020` | 交接时的证据清单 | 历史快照 |

除Python缓存与明确排除的文献原包外，上述近期执行目录保留实际文件，包含失败和首次重试记录。FoldX结果表可以纳入，授权执行程序不纳入。

## 2. 历史与来源审计的精选数据

| 目录 | 纳入范围 |
|---|---|
| `feedback_control_1ycr_20261006` | 根目录实验记录/候选序列/协议与`analysis`直接汇总文件；不含完整结构树与GraphPep输入复制树 |
| `affinity_diagnostic_mdm2_20261006` | 原始亲和力表、冻结/修订记录及`analysis`逐模型、逐序列、方向和相关性表 |
| `foldx_affinity_mdm2_20261006` | 协议、工具清单、输入与执行核对及`analysis`逐模型/方向表；不含整套RepairPDB/AnalyseComplex原始目录 |
| `backbone_diagnostic_20261007` | 作业/源清单、执行记录及`analysis`逐结构/逐配对/能量项/几何/seed表和图；不含240结构及完整FoldX原始树 |
| `independent_affinity_validation_20261007` | 根目录候选、来源审计、原始文献Table2数值和配对记录；保留来源缺项，不重称独立成功 |
| `independent_direction_audit_20261020` | 根目录PDZ/Homer1表，followup、Homer1筛查和PDZ重叠元数据；不纳入论文下载、附件整包及全文提取文本 |
| `bh3_bclw_source_audit_20261020` | 根目录Kd表、输入表示、误差审计、来源与重叠判断；不含raw文献目录 |
| `bclw_second_background_audit_20261020` | 根目录来源测量与背景审计；不含下载/全文提取子目录 |
| `bh3_second_background_audit_20261020` | 根目录后续审计结论；不含文献原包 |

这部分是明确的精选快照。报告链接和原manifest可能指向未纳入Git的历史原始文件；完整数值再审计仍需原工作区。较早pilot、全池生成、完整匹配诊断等目录默认留在本地。

## 3. 未纳入Git的内容

- 任何模型权重与环境目录、授权FoldX程序及其安装文件。
- Python字节码和缓存、传输用重复Base64脚本、临时文件、会话启动提示。
- 论文PDF/XML/Office原包、归档压缩文件和原始全文提取目录；来源URL、哈希和派生数值审计保留。
- 大部分早期AFM结构树、GraphPep特征/输入复制目录及240结构原始树。
- 上游参考notebook没有作为本项目自有代码纳入。

`.gitignore`仅控制Git收录，不删除本地文件。没有为节省仓库体积重写原始结果、清单或冻结配置。

## 4. 如何理解并核验清单

`repository_data_manifest.json`是此次Git数据边界的逐文件清单；每条包含相对路径、字节数和SHA256。README提供标准库核验命令。

原研究目录中的`artifact_manifest.json`、`run_result_manifest_v1.json`、`preservation_before.json`等原样保留。部分原清单含被排除的缓存、文献原包或绝对路径，因此克隆仓库不保证能直接通过每一种旧本地审计入口。不能删除旧清单中的缺项来伪造完整恢复。

本次使用`.gitattributes`禁用换行转换，并核对Git暂存对象与原文件字节，使已经纳入的冻结文件和结果可以原样传递。克隆后先核对本次数据清单，再根据具体研究任务决定是否恢复未纳入资产。

## 5. 来源与标签

数据来自项目执行、公开论文测量、UniProt、RCSB PDB和SKEMPI等。具体来源、测量对象、输入化学和缺项见各审计表。第三方材料依原许可使用，数据库重复条目、反向边和不同受体不能当作新增独立样本。

这份数据快照不含本项目新设计肽的实验Kd；也没有因此建立完整独立benchmark。H原型和Bcl-w结构诊断保持各自的未验证状态。
