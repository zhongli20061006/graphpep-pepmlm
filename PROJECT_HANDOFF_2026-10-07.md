# Peptide design 闭环研究：工作交接

更新时间：2026-10-07（北京时间）。本文件用于新会话接续工作；数据结论以对应实验的冻结协议、原始表与审计结果为准。

## 1. 接手时首先知道的事

**工程流程已跑通，但“反馈引导能稳定提高真实亲和力”的科学命题尚未验证。** 当前最重要的瓶颈是优化目标与选位规则的有效性，继续扩大生成数量无法替代这项验证。

研究目标：在低算力、小规模、不训练大型模型的条件下，检验下面的流程是否有用：

```text
受体序列 + 父代肽
    → PepMLM 局部 mask 重采样产生子代序列
    → AFM / ColabFold 预测每条子代与受体的复合物
    → 结构质量筛选
    → GraphPep 复合物评分及逐残基接触反馈
    → 与父代及匹配随机策略比较
    → 有效反馈用于下一轮局部重生成
```

已经完成一轮“反馈选位 → 重生成 → 子代重新折叠 → 评分”的工程验证。**没有完成子代成为新父代的多代迭代，也没有新设计肽的实验 Kd。** 后续构想是保留 GraphPep 的结构/接触评价用途，另行验证亲和力优化目标，再恢复生成闭环。

最初靶标 NCAM1 的结合位点未知，生成了 1000 条 10-aa 候选。小规模结构试验质量不足后，转向有 PDB 复合物的 MDM2–肽体系进行方法诊断。两个体系不能合并解释成 NCAM1 已验证成功。

## 2. 工作方式与约束

- 使用中文解释，清楚区分“实际完成”“计划”“计算代理改善”和“实验验证”。用户希望直接推进必要工作，不要频繁询问已授权的常规操作。
- RTX 5060 Laptop GPU，约 8 GiB 显存；Windows + Ubuntu WSL。用户描述本机 CUDA 12.8；实际包版本以环境检查和运行记录为准。
- 不需要训练大型模型，不默认购买云算力，不默认折叠全部 1000 条。
- GraphPep 已部署并验证，**不要再次下载或重装本体**。用户指定位置是 `\\wsl.localhost\Ubuntu\home\g\GraphPep_v1.3`。
- 用户已从官网下载授权 Linux FoldX 包，解压后的现有 executable 可以使用；不要重复要求注册下载。
- 保留旧实验、原始结构、协议与审计记录。新实验使用新输出目录，避免覆写冻结结果。
- 旧会话曾做过代码/科学审查，但新会话不应假设旧代理仍在运行，或依赖它们的未落盘上下文。

## 3. 目录、环境与运行入口

|用途|位置|
|---|---|
|Windows 主工作区|`E:\Users\graphpep与pepmlm`|
|同一工作区的 WSL 路径|`/mnt/e/Users/graphpep与pepmlm`|
|早期 WSL 项目与模型目录|`/home/g/pepmlm_pipeline`|
|GraphPep 仓库|`/home/g/GraphPep_v1.3`|
|PepMLM 权重|`/home/g/pepmlm_pipeline/models/PepMLM-650M`|
|GraphPep / PepMLM Python|`/home/g/miniconda3/envs/graphpep_env/bin/python`|
|ColabFold Python|`/home/g/miniconda3/envs/colabfold_exp0/bin/python`|
|ColabFold CLI|`/home/g/miniconda3/envs/colabfold_exp0/bin/colabfold_batch`|
|FoldX executable|`/mnt/e/Users/graphpep与pepmlm/tools/foldx/foldx5_Linux/foldx_20261231`|
|FoldX 原包|`tools/foldx/foldx5_Linux_0.zip`（相对主工作区）|
|方法计划|`docs/superpowers/plans/`|
|主要源码 / 运行脚本 / 测试|`pipeline/`、`scripts/`、`tests/`|

上述 WSL 项目、GraphPep、两个 Python、ColabFold CLI、PepMLM 权重目录于交接时只读确认存在。未在本次交接中重新启动 GPU 推理，也未全面重新核验依赖版本。

早期脚本有硬编码 `/home/g/pepmlm_pipeline`，后续部分脚本按文件所在位置定位项目。**两个项目目录不是可随意互换的单一目录。** `scripts/sync_pipeline.sh` 会把 Windows 的 `pipeline/*.py` 复制到 WSL 旧项目并转 LF；这是写操作，不要为了读报告盲目执行。修改后先明确实际运行哪份源码，再做必要同步和版本/哈希记录。

在 WSL 中的轻量检查示例：

```bash
cd /mnt/e/Users/graphpep与pepmlm
/home/g/miniconda3/envs/graphpep_env/bin/python --version
/home/g/miniconda3/envs/colabfold_exp0/bin/colabfold_batch --help
```

Windows 中可通过 `wsl.exe -d Ubuntu -- bash -lc '命令'` 进入。受限环境调用 WSL 可能需要工具级授权；这不代表依赖损坏。PowerShell 读取中文文档请显式使用 `Get-Content -Encoding UTF8`。

已用过的 AFM 配置：ColabFold 1.6.3、`alphafold2_multimer_v3`、通常参数模型 1/2/3、每次一个 seed、3 recycle、`max_seq=256`、`max_extra_seq=512`、无模板、先不 relax，再用统一 CPU Amber 松弛。启用 `--disable-unified-memory` 关闭 unified memory 曾解决本机推理异常；其他主要参数见下列已执行脚本，不能假设适合所有新受体：

- `scripts/run_feedback_control_1ycr.sh`
- `scripts/run_affinity_diagnostic.sh`
- `scripts/run_foldx_affinity.sh`

这些脚本是**历史运行入口**，不是接手后应立即执行的待办清单。ColabFold 的 MSA 检索可能使用共享服务器，本地 AFM 推理不等于完全离线。

## 4. 已完成实验与核心结果

### 4.1 初始 NCAM1 生成池

原始 CSV 在 `/home/g/pepmlm_pipeline/output/generated_peptides.csv`；当前工作区保留 `output/experiment0/generated_peptides_snapshot.csv`。

本次交接从快照重新核对：1000 行、976 条唯一序列、全部为标准氨基酸组成的 10-mer；PPL 最小 6.464576、中位 19.365826、最大 51.138157。

早期运行日志记录总耗时 826 s，fp16 采样从 154.5 s 降至 5.4 s；PPL 与作者逐条实现相对误差约 3.4e-4。ESM-2 转换逐 logit 最大差约 1.35e-4。上述性能/等价性数字是既有运行记录，本次没有重复计算。

阶段一报告记录：位置熵均值 0.342（归一到 0–1），K 占 23.6%；低 PPL Top-100 的 K 升到 31.3%。这是组成与多样性问题，**不等于已经证明局部重采样能改善结合**。PPL 是生成模型的序列指标，不是亲和力。

### 4.2 Experiment 0：本地预测与 GraphPep 输入验证

报告：`output/experiment0/local_trial/local_trial_report_zh.md`；扩展报告：`output/experiment0/expanded_validation/expanded_validation_report_zh.md`。

- 受体：探索性 NCAM1 20–718，699 aa；全长 858 aa。两条候选：PEP00815 `TKKRQAARGK`；PEP00356 `TIQKKLLGIP`。
- 两条候选各 5 个 AFM 参数模型；另做 1YCR / MDM2–p53 的 5 个预测模型及 1 个实验结构正对照。15 个预测模型完成，16 份结构可由 GraphPep 处理。
- NCAM1 两组平均 ipTM 约 0.370 / 0.336，肽 pLDDT 20.76 / 17.76，界面 PAE 19.12 / 20.30 Å；各 5/5 有 <1.5 Å 碰撞，接触集合跨模型不稳定。**不能用于闭环反馈，也不能由此证明实际不结合。**
- 1YCR 正对照平均 ipTM 0.840、肽 pLDDT 90.08、受体对齐后肽 CA RMSD 1.68 Å，恢复约 71.3% 原生 5 Å 接触。
- 正对照预测 GraphPep 分数 −3.495 ± 0.037，实验结构 −3.829；整体评分较稳定，但预测与晶体逐残基反馈 Pearson 为 −0.27 至 0.21，不能认为反馈已经普遍可靠。
- 五个参数模型共用 seed 20261002，**不是五个独立随机 seed**。生成与预测使用的 NCAM1 受体范围还不同。

NCAM1 路线当前暂缓。重新启动需要结合区域或可靠阳性配体的独立依据，并统一生成/预测构建体。

### 4.3 MDM2：反馈重生成试验及复验

父代主要为 1YCR 的 13-mer `ETFSDLWKLLPEN`，受体 MDM2 25–109。后来扩展至 3EQS 的 PMI 12-mer 和 3G03 的 p53-like 12-mer。

|实验目录（均在 output 下）|设计与主要发现|
|---|---|
|`feedback_pilot_mdm2_20261005`|实验结构反馈选位 3 条 vs 匹配随机 3 条；初始代理评分达标 2/3 vs 0/3。只是小样本探索。|
|`feedback_seed_retest_mdm2_20261005`|统一 Amber 后旧 seed 3/3 vs 0/3，新 seed 2/3 vs 1/3；一条引导序列跨 seed 重复，且采样与折叠 seed 同时改变，不算独立新序列复现。不能与未松弛 pilot 直接合并计数。|
|`feedback_panel_mdm2_20261005`|3 个父代 × 4 个策略臂 × 每臂 3 条子代；**36 条子代 + 3 条父代 = 39 条序列，共 117 个预测模型**。实验结构引导 2/9 vs 随机 2/9；预测结构引导 2/9 vs 随机 1/9，未显示可靠普遍优势。|
|`feedback_control_1ycr_20261006`|严格位点池匹配，预测反馈引导 10 unique vs 随机 10 unique。60 个子代模型均过结构门槛，两组代理改善均 2/10；平均配对中位增益 +0.01714 vs +0.02084。没有证明引导更优。|

严格对照保护 F3/W7/L10，每条替换两个非锚点，共同可选池为 `[1,2,4,6,9,11,12,13]`。引导只选 E1/N13，随机覆盖多个位点对；十条引导子代不能当作十次独立选位验证。两组 top-k=4，PPL 仅记录、不筛选。

3G03 中 S12P 曾出现在达标子代中，但该末端残基在实验坐标中未解析。严格 1YCR 对照中达标引导子代含 N13P，达标随机子代含 E12P；这些观察没有实验 Kd 支持，不能声称末端 Pro 普遍增强亲和力。

### 4.4 GraphPep：已知亲和力诊断

报告：`output/affinity_diagnostic_mdm2_20261006/report_zh.md`。

绕过 PepMLM，用 Li et al. 2010 同一 SPR 竞争实验、同一 synthetic MDM2(25–109) 构建体的 10 条已知 12-mer：PMI 与 p53 两个家族各 5 条。两 seed（20261011/20261012）× 三参数模型，共 60 个结构，均通过当前质量门槛。

- 总体 GraphPep vs log10(Kd) Spearman 0.736 / 0.648。
- PMI 内 1.000 / 1.000；p53 内 0.100 / 0.000。总体排名不能替代局部突变排序。
- 七项 ≥3 倍的强效突变方向，两 seed 均 5/7 正确；P27A、L22A 都反向。
- 预设继续门槛失败；旧规则漏掉无接触 N8A 和高接触置信度 P27A 的已知收益。
- PMI WT / N8A 的实验结构共同范围对照也不能稳定复现亲和力方向。它只有两个序列、且存在末端裁剪，不能据此完成各模块的严格因果归因。

此结果支持**评分/接触置信度与亲和力优化目标存在错配**。不是 PepMLM 部署错误的直接证据，也不能只归咎于 AFM。GraphPep 可继续作为结构/接触评价工具，但尚不能当作已经校准的跨序列亲和力目标。

### 4.5 FoldX：独立代理目标开发筛查

报告：`output/foldx_affinity_mdm2_20261006/analysis/report_zh.md`；判定：同目录 `decision.json`。

使用同一已知十肽面板的 60 个 AFM 模型，FoldX 评分全部成功、独立输入与原始输出审计通过。

|seed|总体 / PMI / p53 Spearman|强效方向|预设门槛|
|---|---|---|---|
|20261011|0.976 / 1.000 / 0.900|6/7，P27A 错误|FAIL|
|20261012|0.867 / 1.000 / 1.000|7/7|PASS|

两个 seed 联合 `continuation_pass=false`。P27A 配对增益约 −0.114 / +0.046 kcal/mol，方向依赖 seed。整体排名不错，仍不满足稳定局部优化要求。这里用的是 `AnalyseComplex` 的 A/B Interaction Energy，**不是已验证的 ΔG(binding)、实验 Kd 或精确 ΔΔG**。

### 4.6 固定骨架与释放肽主链诊断

固定骨架报告：`output/foldx_matched_diagnostic_20261006/analysis/interpretation_zh.md`。24 个任务、120 个 WT/突变配对、240 次 AnalyseComplex；P27A 30/30 配对反向。它是结合态固定骨架敏感性检查。

最新主链诊断报告：`output/backbone_diagnostic_20261007/analysis/interpretation_zh.md`。

- WT/P27A/L22A，同一起点的 fixed vs flexible，对齐两 seed、三参数模型、五重复。240 个结构有效，120 个条件内配对，60 个 flexible−fixed 配对差；48 个试跑结构**包含在 240 个内**。
- 独立审计复核 480 张原始表、8160 个数值字段及输入/环境/代码/产物哈希；技术与数据完整性通过。
- P27A：fixed/flexible 下 seed 20261011 均错误方向，20261012 均不确定；释放主链没有恢复稳定预期收益。L22A 各条件方向正确。
- P27A 的 flexible−fixed 配对效果 +0.1044 / +0.0917 kcal/mol；是先对相同重复作差后聚合，不能用表中汇总中位数直接相减替代。
- 肽位置 8–11 的 alpha-like 残基中位数 1→2，但 WT 也 1→2，突变特异配对变化中位数为 0；O7···N11 6.20→6.07 Å，没有结构达到 3.5 Å 几何代理阈值。
- 严格可动原子向量力 RMS ≤10 kJ/mol/nm：35/240 通过、205/240 未通过；内部优化器结束不等于此严格验收通过。
- 原有独立 AFM P27A 结构已有末端变化，受体对齐后 B8–12 CA RMSD 约 2.7–2.9 Å。不能说 AFM 完全没有构象响应。
- OpenMM 8.6.1 CPU、ff14SB/OBC2、局部最小化；没有游离肽/受体系综，也没有完整热采样。这不是 MD 自由能验证。当前仍 `continuation_pass=false`。

这轮完成后没有新研究计算。P27A 仍是有用反例，但目前不宜继续无界扩大相同局部最小化。

## 5. 必须沿用的指标口径与解释边界

1. GraphPep 原始分数越低（通常越负）越有利。增益 = 父代 − 子代，正数表示代理评分改善；不能叫作实验亲和力提升。
2. 严格 10v10 主指标是 `median(parent_model_score − child_model_score)`，使用相同参数模型、质量通过的配对。`median(parent) − median(child)` 是不同公式；旧 panel 用过后者。严格对照的事后替代公式会变成 3/10 vs 2/10，但不得替换预设 2/10 vs 2/10 结论。
3. 常用结构门槛：ipTM≥0.70、肽 pLDDT≥70、界面 PAE≤10 Å、有接触、链间最小重原子距离≥1.5 Å；实际 PAE 聚合口径按各 protocol。至少 2 个相同参数模型父子配对通过。
4. 操作性子代改善：配对中位增益 >0.05，且至少两个配对方向改善。0.05 没有实验能量单位，也不等于统计显著性。
5. 接触置信度不是残基能量贡献，也不是突变收益梯度；无接触不等于“不值得改”。缺失反馈与数值为零分开记录。
6. AFM 参数模型/随机 seed、晶体副本、FoldX 数值重复都不是生物学重复。Amber 后不能重新解释原 AFM pLDDT/ipTM 为已更新的置信度。
7. 实验亲和力面板是 12-mer p53 `ETFSDLWKLLPE`，重生成父代是 13-mer `ETFSDLWKLLPEN`。P27A 对应 p53 编号 27、12-mer 肽位置 11；不要混淆构建体与编号。
8. 已知十肽已经用于方法开发，不能再称为独立盲测。训练集/结构重叠尚未审计；独立数据验证也不自动等于模型训练集无重叠。

GraphPep 输入适配由 `pipeline/step3_graphpep.py` 处理：各序列对应各自真实复合物；先核对链序列并确定受体/肽角色，保留正确坐标与编号。不要将不同序列共享一个 multi-MODEL 拓扑或 ESM 特征，不要用确定性摆放替代复合物预测。此前兼容补丁位于项目适配层，原 GraphPep 核心源码哈希在审计中保持不变。

## 6. 证据索引与需警惕的旧文件

以下路径除标明外，均相对 `E:\Users\graphpep与pepmlm`。读取顺序建议先最新解释、再相应协议/审计/明细。

|问题|优先证据|
|---|---|
|最新结论和主链限制|`output/backbone_diagnostic_20261007/analysis/interpretation_zh.md`、`protocol.json`、`run_record.json`、`analysis/independent_audit.json`|
|FoldX 亲和力开发筛查|`output/foldx_affinity_mdm2_20261006/analysis/report_zh.md`、`analysis/decision.json`、`analysis/independent_score_verification.json`|
|GraphPep 的亲和力/选位反例|`output/affinity_diagnostic_mdm2_20261006/report_zh.md`、`analysis/mutation_directions.csv`、`analysis/feedback_vs_mutation.csv`、`analysis/independent_verification.json`|
|严格 10v10 对照|`output/feedback_control_1ycr_20261006/report_zh.md`、`protocol.json`、`experiment_record.csv`、`analysis/aggregation_sensitivity.csv`|
|多父代对照|`output/feedback_panel_mdm2_20261005/report_zh.md`、`experiment_record.csv`、`analysis/verification.json`|
|初始结构验证|`output/experiment0/expanded_validation/expanded_validation_report_zh.md`、`validation_summary.json`（后者位于 expanded_validation 下）|
|生成池快照|`output/experiment0/generated_peptides_snapshot.csv`|
|汇报 Word|`output/research_report_zh/PeptideDesign_ResearchUpdate_2026-10-07.docx`|

科研汇报 Word 已交付桌面：`C:\Users\g\Desktop\PeptideDesign_ResearchUpdate_2026-10-07.docx`。结构为约 11 页，包含方法、当前瓶颈和后续建议；实际 Word 分页未经原生渲染核验。本次交接查到工作区副本存在，319604 bytes。

**旧状态文件不能覆盖较新报告：**

- `/home/g/pepmlm_pipeline/output/stage1_report.md` 是 interim，Top-20 数为 0，BindEvaluator 分数未就绪，结构 pending。其“GraphPep higher better”写错，应按后续适配和原始分数采用 lower better。其他旧措辞也不要直接当作科学结论。
- `output/experiment0/README_zh.md` 与 `local_verification.json` 记录的是云端准备阶段“尚未启动”；后来已经完成本地 AFM，不能据此判定当前未折叠。
- `tools/foldx/README.md` 仍描述待用户下载的准备阶段，实际授权 executable、60 模型评分和后续诊断已经完成。
- 个别旧实验报告写“亲和力诊断尚未执行”，仅代表其生成时点，后续实验已执行；每个报告的阴性结论和协议仍应保留。
- BindEvaluator 初始下载/打分没有在这次接手核查中找到最终 Top-20 交付证据，不应当作已经完成，也不建议把补齐它设为当前默认主线。
- 本次没有核验所有后台任务或运行中的下载进程；不能仅凭历史 ETA 推断当前仍在下载。

## 7. 下一阶段建议：先校准目标，再扩展闭环

下面是建议方向，**尚未执行**。新会话应先阅读与核对，再把方案落实为可评审的协议/数据表。

### 优先任务 A：独立亲和力验证方案与小型数据集

1. 整理已用十肽为 development set，冻结当前方法、结构质量处理、失败/无接触处理、汇总公式和预期方向规则。
2. 从关键原始论文和结构数据库寻找未参与当前选法的、同一受体/同一 assay 内可比较的肽及局部突变数据；优先有实验复合物、构建体/末端化学明确、已知增强和减弱突变的体系。MDM2 可作一部分，但不要只选符合当前评分预期的突变。
3. 建立候选表：target/构建体、肽序列及修饰/长度、测量方法、Kd及误差/删失、来源表格、结构来源、parent–mutant关系、是否用于开发、训练集重叠核验状态。
4. 在看新分数前冻结判据。重点看**同家族排序、局部改善方向、两个折叠 seed 的稳定性**，而不是仅 pooled rho。小样本描述性结果与置信/显著性分开。
5. 优先测试现有 FoldX 等低预算代理的独立有效性；明确没有游离态构象采样的局限。不要默认简单 MM/GBSA 或局部松弛已经解决熵/螺旋问题。
6. 报告方法开发 holdout 与真正外部盲测的区别；新集合一旦用于调参就不再是未触碰验证集。先完成数据审计和预注册式协议，再决定是否值得新预测。

### 条件任务 B：重新做少量设计与匹配随机对照

仅当目标在独立数据上有足够支持后，选 1–2 个结构与结合机制明确的体系，每臂 10–20 条 unique 子代；统一可编辑位点、锚点保护、突变数、采样预算、受体/MSA、折叠及松弛条件。采样 seed 与折叠 seed 分别记录/控制。允许无接触位点进入预先声明的候选池，检验选位而非只检验一个末端位点对。

GraphPep 负责结构/接触信号，独立亲和力目标负责优化终点；设匹配随机、必要的序列模型基线和负对照。只有一轮有稳健证据后，再让子代成为父代做第二代。

### 条件任务 C：实验验证

准备少量“父代 + 引导子代 + 随机子代 + 负对照”的同 assay、同构建体测量面板，主结果以实验 log(Kd) 或明确亲和力变化衡量。当前没有合成、SPR/ITC 订单或新测量授权，先产出可审阅清单与设计，不联系第三方或付费。

### NCAM1 的重启条件

独立结合位点/已验证配体/目标相关阳性对照，统一受体构建体，少量多 seed 结构达到界面可信度要求。当前不恢复全池折叠。

## 8. 新会话的首轮交付建议

先简要确认项目到哪一步；读本文件和最近三份诊断报告。只读核验环境、协议与结果存在，避免自动重跑历史 shell 脚本。随后产出一份“独立亲和力验证协议草案 + 数据来源审计表”，列出 1–2 个可行体系、计算预算、明确主指标、进入/停止条件与缺失证据。对原论文未提供的数值留空或标为待核验，不推测。

当前用户最新授权是写交接文档与提示词，并未要求本会话立即执行下一轮研究实验。新会话收到继续指令后再推进对应工作。

## 9. 关键文献与原始结构

- GraphPep：Tao et al., 2025，An interaction-derived graph learning framework for scoring protein–peptide complexes，DOI 10.1038/s42256-025-01136-1。https://www.nature.com/articles/s42256-025-01136-1
- 实测十肽面板：Li et al., 2010，Systematic Mutational Analysis of Peptide Inhibition of the p53-MDM2/MDMX Interactions，DOI 10.1016/j.jmb.2010.03.005。https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/
- p53 构象机制参考：Higbee et al., 2024，Structural Adaptation of Secondary p53 Binding Sites on MDM2 and MDMX，DOI 10.1016/j.jmb.2024.168626。https://pmc.ncbi.nlm.nih.gov/articles/PMC12624572/ （构建体与当前 12/13-mer 不同，不能直接套用。）
- ColabFold：Mirdita et al., 2022，DOI 10.1038/s41592-022-01488-1。https://www.nature.com/articles/s41592-022-01488-1
- 原始结构：https://www.rcsb.org/structure/1YCR 、https://www.rcsb.org/structure/3EQS 、https://www.rcsb.org/structure/3G03 、https://www.rcsb.org/structure/3LNZ

新增方法与数据来源需在新会话按原始论文、官方文档重新查证；本索引不是最新数据库的全面检索结果。
