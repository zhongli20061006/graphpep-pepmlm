# a1-syntrophin PDZ 暴露与重叠审计（2026-10-20）

本轮只读取原有文件，并把公共 RCSB 序列/元数据和新审计记录写入本轮目录。未运行折叠、评分、模型训练、安装或权重下载。旧目录不修改。

## 结论与边界

原论文把三个母肽用于自己的训练数据，并以单点变体检验其自己的预测器；这不自动把它们变成本项目或 FoldX/GraphPep 的开发数据。10 月 8 日本项目已有候选来源/选法讨论，因而本轮不是从未见过候选身份的前瞻盲选；本轮读取标签、未看任何新评分，后续若使用，应称资料审计后的外部评价，披露选体系过程。

GraphPep 所附训练/验证/测试列表没有 2PDZ exact ID；2509 套原生受体/肽结构中也没有三个母肽 core10 或 NNG+core 的 full13 exact/substring 命中。2PDZ/tutorial 的受体参考序列为 **86 aa**，不是此前估计的 87 aa；本地也没有该 86 aa 全长 exact receptor sequence 命中。它是参考/教程序列；Supplementary Table 2直接声明的assay construct为UPI0000566AF6的90–164（75 aa）。已从UniParc原始504 aa序列恢复这75 aa，为86 aa参考序列的第12–86位；两者角色不同，本身不是文献矛盾。75 aa为`GGLGISIKGGRENKMPILISKIFKGLAADQTEALFVGDAILSVNGEDLSSATHDEAVQALKKTGKEVVLEVKYMK`。额外扫描2509个本地受体亦无75 aa exact/substring命中；N端Trx/His6及接头完整序列尚未恢复，因此完整表达构建体输入仍未解决，不能用86 aa代替75 aa+tags。

有明确的 **PDZ 同家族数据暴露**：用该 86 aa 在 RCSB 蛋白序列服务检索（identity >=0.30、E-value <=0.01），得到 265 个 polymer entity、263 个 PDB entry；与 GraphPep 所附 Training_set 交叉为 39 个 target，其中 train 35、valid 4。最高相似度是 50.0% sequence identity，query 覆盖 84/86=97.67%，例为训练集 4E34、4JOE 等 CAL PDZ 与验证集 4JOP。另有 Scribble PDZ1、PSD95 PDZ3 等。**这些是同折叠/同家族受体，不是 alpha1-syntrophin/SNTA1 同受体泄漏的证据**。39 条不是39个独立受体，更不等于39份亲和力突变标签。 对这39个交叉polymer entity逐条读取RCSB蛋白名/gene元数据，39/39成功、无SNTA1或alpha1-syntrophin文本命中；完整来源保留在family_name_gene_audit.json。

未证明所附 splits 就是当前待用 GraphPep checkpoint 的全部实际训练来源。ESM-2 与 AFM/AF2 的确切预训练成员未审计，不能称预训练独立。可公开披露家族暴露；不能以 exact ID 未命中声称完全无重叠。

本轮停止 pilot 的主要原因是原始表中没有达到冻结 **>=2 倍增强**阈值的变体，不是 PDZ 同家族在 GraphPep 训练中出现；完整构建体问题亦尚需核清。

## 本地读取范围

- Training_set/list：2231 target；Training_set/split_data/train.csv：2007 target / 32666 data rows；valid.csv：224 target / 3635 data rows。根目录 split_data 不存在，真实路径在 Training_set 下。
- 测试 manifest：ADCP99 99；LEADS-PEP 53；AF3 list 64 / list_all 121；Local62 62。
- 原生序列扫描：Training_set 2231、ADCP99 99、LEADS53、AF3 64、Local62 62，共2509套 / 5018个受体与肽输入文件。只读第一 MODEL 的 CA 记录，按 chain/residue ID去重；MSE映射M；Local62受体使用 recb.pdb。
- 未扫描巨型 decoy 集、逐模型 result JSON或模型权重。追加核对已覆盖pdz_ligand_input_audit.csv的全部18行（3母肽+15变体），core10及NNG+core13共36个唯一查询均无exact或query-in-native命中；同时做native片段（>=6aa）in-query补查亦无命中。查询表源SHA256已记录于all18_pdz_peptide_overlap.json，化学修饰与不完整建模残基仍不受此结果排除。
- 限定项目报告、Markdown、manifest/audit JSON、CSV/TSV/FASTA等1913文件。仅 output/reward_direction_benchmark_20261008/protocol_draft_zh.md 与 candidate_source_audit.csv命中候选词；未找到更早日期输入/报告命中，也未找到候选新评分产物。没有git历史，目录日期/当前文件系统不构成完整时间溯源证据。
- E:/AGENTS.md 是FastAPI专用说明，已读取，本审计无FastAPI实现；工作区、output与WSL GraphPep祖先未发现其他适用AGENTS。

## 可核验的证据

[pdz_overlap_audit.json](pdz_overlap_audit.json) 汇总结论、列表计数和SHA256。pdz_overlap_sources/structure_input_hashes.json列出5018个实际读取结构的路径、大小、SHA256；graphpep_sequence_scan.json记录核心/全长肽查询、未命中和初筛候选；rcsb_sequence_search_30pct_verbose.json保留服务器给出的逐命中identity/覆盖字段；split_family_matches.json与sequence_identity_summary.json可复算train/valid交叉。workspace_match_audit.json包含本地命中文件和行号。审计JSON的sources列记录来源文件哈希。

RCSB full_text搜索“syntrophin”还曾返回测试集1EG4；逐实体核对其实际为dystrophin WW–beta-dystroglycan，属于全文检索的假阳性，未将其记为syntrophin受体重叠。

原始来源：[Chen et al. 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2655215/)、[2PDZ](https://www.rcsb.org/structure/2PDZ)、[2PDZ FASTA](https://www.rcsb.org/fasta/entry/2PDZ/display)、[4E34](https://www.rcsb.org/structure/4E34)、[6MS1](https://www.rcsb.org/structure/6MS1)、[1BE9](https://www.rcsb.org/structure/1BE9)。2PDZ database description可标记GVKESLV，但canonical polymer FASTA返回可建模的KESLV；两者都不是本实验的TAMRA-NNG+10mer，应保留该差异。


## 窄范围补查：Homer1 EVH1（与PDZ来源分开）

从本轮homer1_screen/homer1_single_mutants.json的11个单点配对提取15个唯一parent/variant肽插入序列，查询同一2509个native peptide文件；exact、query-in-native、native片段>=6aa-in-query均零命中。查询源路径、SHA256、每条序列与结果保留于pdz_overlap_sources/homer1_narrow_overlap.json。**未恢复完整生物素化SUMO融合构建体，不宣称完整assay sequence无重叠。**

参考结构1DDV exact ID不在本地Training_set/train/valid及其他已查列表，却在GraphPep所附Local_62_docking_decoys/list中命中1DDV。这是分发测试集结构曝光，不等于亲和力突变标签训练曝光；本窄补查没有下载更多家族结构、没有做Homer1同受体序列或家族重叠审计。

限定的此前工作区md、manifest/audit JSON、CSV/TSV/FASTA/FA搜索未见Homer1、1DDV或PMC11237882命中（排除本轮目录）；这不证明机器/预训练全集未暴露。本轮已见该来源及11对标签，未看新评分，不称从未见标签的双盲。本补查不推断原论文的模型开发用途。

原始来源：[Homer1 EVH1论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11237882/)。
