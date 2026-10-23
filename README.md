# 低算力蛋白—肽设计闭环：PepMLM / GraphPep

这是一个研究工程工作区，保存生成、复合物预测、候选评分、逐残基反馈和局部重生成的实现，以及支持或否定各项假设的阶段证据。

**当前状态（2026-10-23，Asia/Shanghai）：工程闭环已跑通；尚未证明反馈优于匹配随机，也未证明新设计肽的实验亲和力改善。H原型未观察到额外开发收益；Bcl-w的小型MSA输入审计已通过，亲本双seed结构诊断已准备、未启动。**

本仓库保存代码和选择纳入的本地阶段数据。模型权重、软件环境、授权FoldX程序、论文下载整包和大部分历史结构树仍由原工作区保存。数据边界和字节校验见[数据目录](docs/DATA_CATALOG.md)及[数据文件清单](docs/repository_data_manifest.json)。

本次选择纳入808个数据文件，共34,715,464字节（约33.1 MiB），采用普通Git保存，未启用Git LFS。

## 1. 目标与当前研究问题

长期目标是在约8 GB显存条件下，进行小规模、可证伪的肽设计验证：

```text
PepMLM生成
    ↓
ColabFold / AFM复合物预测
    ↓
结构质量、候选评分及逐残基反馈
    ↓
局部mask重生成
    ↓
重新预测评价，并与匹配随机对照比较
```

不训练大型模型。当前首先要解决的是：什么信号能够可靠判断同一受体、同一assay、同一亲本肽单点突变的亲和力方向。只有候选reward获得支持后，才继续投入反馈设计和多代优化。

没有完成的目标包括：新设计肽实验Kd、完整ΔG(binding)、多代优化，以及NCAM1泛化。

## 2. 已完成的关键工作与证据

| 阶段 | 已完成 | 结果及解释 | 记录入口 |
|---|---|---|---|
| 严格反馈对照 | 引导10条与匹配随机10条；子代预测评价完成 | GraphPep代理改善各2/10，未显示引导优势 | [严格10v10](output/feedback_control_1ycr_20261006/report_zh.md) |
| GraphPep方向诊断 | 已知亲和力MDM2开发面板 | 存在局部方向反例；总体相关性不能替代亲本内方向 | [方向诊断](output/affinity_diagnostic_mdm2_20261006/report_zh.md) |
| FoldX诊断 | 同一面板跨seed评价 | 总体排序较好，但P27A跨seed不稳，联合continuation_pass=false | [FoldX结果](output/foldx_affinity_mdm2_20261006/analysis/report_zh.md) |
| 主链敏感性 | 240个结构诊断已完成 | 35/240通过严格力RMS，未稳定恢复P27A；不重复此批最小化 | [主链解释](output/backbone_diagnostic_20261007/analysis/interpretation_zh.md) |
| Bcl-w六结构pilot | 三条肽×两个固定seed；预测、放松和评分完成 | 0/6通过冻结结构门，只作工程失败诊断 | [pilot结果](output/bh3_bclw_six_complex_pilot_20261020/run_result_report_zh.md) |
| Bcl-w kernel诊断 | 亲本×两个seed，关闭fast kernels | 0/2通过全部诊断门；不能把原因单独归给fast kernels | [kernel结果](output/bclw_parent_kernel_ablation_v1/report_zh.md) |
| H技术原型 | 60条逻辑模型输入完成；没有折叠新子代 | 数值检查通过，建议会变化；未验证亲和力 | [原型结果](output/receptor_context_mutation_probe_v0_20261022/report_zh.md) |
| H开发比较 | 七对PMI点效应方向比较；19条新逻辑输入 | d和H均7/7，没有观察到额外收益；不是正式独立方向准确率 | [PMI比较](output/receptor_context_pmi_development_comparison_20261022/report_zh.md) |
| H分歧审计 | 复算已有361项替换读数 | 47项严格反向、37项H拒绝；84项分歧无完整匹配实验Kd配对 | [分歧审计](output/receptor_context_disagreement_audit_20261022/report_zh.md) |
| 最新MSA审计 | 93条已有候选完成筛选与A3M检查 | 22条非query同源序列、4簇、Neff≈4.01、最低列覆盖82.61%；输入门通过 | [MSA报告](output/bclw_msa_input_audit_20261023/report_zh.md) |

低置信度不证明实验不结合。输入门通过也不证明预测结构会恢复。模型概率、几何接触和结构能量均不等于实验Kd或完整结合自由能。

## 3. H原型具体做什么

H是使用现有PepMLM的推理时对照方法，没有新训练。它尝试扣除对受体顺序不敏感的替换偏好，直接比较局部位置与替换氨基酸。

对受体R、亲本肽P的位置i，只mask该位置：

```text
d_target = log p(替换残基 | R, P_mask_i)
         - log p(原残基   | R, P_mask_i)

H = d_target - (d_shuffle1 + d_shuffle2) / 2
H1 = d_target - d_shuffle1
H2 = d_target - d_shuffle2
```

两个shuffle受体与真实受体长度、组成相同，序列顺序按固定种子打乱。它们是人工上下文扰动，不是实验证实不结合的受体，可能偏离模型训练分布。

- 开发方向比较：d以±1e-4 nats为数值容差；H1/H2均正预测增强、均负预测减弱，否则拒绝。拒绝保留分母。
- 候选建议阶段另要求d_target、H1、H2都大于容差，随后按H排序，每个位点最多一个替换、每亲本最多六条。
- 1e-4只是数值容差，不是生物学效应阈值。
- 有完整输入匹配的本地实验配对共10项，H/d全部同向；其中七项进入了有条件开发比较。
- 当前结论是“额外收益未证实”，保留本版原型、暂停扩算。不能修改控制或阈值后，再用同一面板宣布成功。

原始概率分布、复用来源、冻结配置、方向表、拒绝及失败记录都在三个`receptor_context_*`目录中。

## 4. 冻结评价口径

正式规则的唯一入口为[frozen_direction_criteria_v2.json](output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json)。文件中的旧状态按生成时点解释；已有执行状态以上表和最新原始记录为准。

- 实验方向使用`Kd_variant / Kd_parent`：≤0.5且95%不确定区间全低于1为增强；≥2且区间全高于1为减弱。
- 需要明确误差种类、重复次数及重复层级。检测不到结合时保留检测限，不编造点Kd。
- 不跨受体、assay、亲本背景拼接配对；Kd、Ki、IC50不能直接混作同一种标签。
- 正式单体系面板要求至少10对、每类至少4对、至少两个亲本背景，另有BA、类别召回和seed稳定性门。这是本项目冻结门，不是所有探索任务的通用样本量规定。
- 已正式合格的实验标签保留失败、无接触、低质量、tie及算法拒绝，不能删去以提高准确率。
- MDM2/PMI/p53已暴露数据属于开发或压力测试；不能称盲测。完整独立benchmark仍未建立。

候选评分方向如下：

| 方法 | 正向差值定义 | 平局/解释 |
|---|---|---|
| GraphPep | parent − variant | 越低越有利；±0.05 score unit平局带；正值只表示代理改善 |
| FoldX | E_parent − E_variant | ±0.05 kcal/mol平局带；不是完整ΔG(binding) |
| ESM2序列PLL | PLL_variant − PLL_parent | 肽序列先验，无受体条件 |
| 接触基线 | contacts_variant − contacts_parent | 几何描述；须明确原子对或残基对 |
| H / d | 上节冻结定义 | nats数值偏好，尚未建立亲和力效力 |

旧v2接触基线写重原子对，而pilot使用唯一残基对；旧PAE协议也有中位数与最大值的差异。新实验必须明确冻结定义，不能把这些列当成同一指标。

## 5. 最新Bcl-w输入与下一步边界

当前核心输入：

```text
受体：已有Q92843全长193 aa中的1–164，133位保持Q
亲本：RPEIWITQELRRIGDEFNAYYAR（23 aa）
MSA：受体query + 22条非query同源序列；肽仅有query
```

22条同源记录中2条reviewed、20条unreviewed。它是小型来源审计比对，不是完整MMseqs搜索，也没有独立证明每条都是严格正交基因。多个同分最优比对已记录，主输入按事前约定取第一个；补充检查没有改变主选择。

模型省略assay的N端荧光团、C端酰胺化及未恢复的受体c-myc/TEV完整接头，仍是核心代理输入。相关BH3结构没有被当成精确assay结构。

已准备协议：[BCLW_PARENT_LOCAL_MSA_DIAGNOSTIC_V1](output/bclw_msa_input_audit_20261023/prepared_diagnostic_protocol.json)。

| 项目 | 冻结内容 |
|---|---|
| 问题 | 这份受体MSA及必要读取/有效深度变化，能否让亲本预测在两seed下恢复工程可用性？ |
| 对照 | 已完成no-fast-kernel亲本分支，复用旧结果，不重跑 |
| 规模 | 亲本×20261008/20261009，最多2个原始结构及2个CPU放松结构 |
| 预算 | GPU标记子进程累计wall≤1800秒；CPU放松≤900秒；输出≤0.5 GiB |
| 原质量门 | ipTM≥0.70；肽均值CA pLDDT≥70；全双向链间PAE最大值≤10 Å；放松后最小重原子距离≥1.5 Å；有接触、序列编号正确 |
| 附加诊断门 | 受体均值CA pLDDT≥70 |
| 结果判定 | 2/2全部通过仅进入单独下一步评审；1/2或0/2为HOLD；执行失败或超预算即停 |
| 当前状态 | 已准备，未启动；不自动续跑突变、评分或闭环 |

已安装ColabFold的`single_sequence`会丢弃同源MSA并限制有效深度为1/1。因此准备的本地A3M方案还须使用非single读取模式。这是干预的一部分，不能描述成只替换文件。协议中的`mmseqs2_uniref_env`是CLI枚举，**不授权远端MSA服务**；本地A3M解析路径已检查，未来启动仍需复核实际输入与保存配置。

## 6. 目录结构与代码入口

```text
README.md                              当前总入口
PROJECT_HANDOFF_2026-10-20.md            历史全貌交接；较晚执行记录优先
PROJECT_HANDOFF_2026-10-07.md            更早历史交接
pipeline/                              生成、预测、反馈和诊断代码
scripts/                               历史运行/环境辅助脚本
tests/                                 原有工程与诊断测试
docs/DATA_CATALOG.md                    本次纳入数据的范围、来源和限制
docs/repository_data_manifest.json      纳入Git的数据逐文件字节数与SHA256
output/receptor_context_*/              H原型、开发比较与分歧证据
output/bclw_msa_input_audit_20261023/    最新输入审计、A3M及未来方案
output/bh3_bclw_six_complex_pilot_20261020/
output/bclw_parent_kernel_ablation_v1/  已完成且失败的结构试验记录
```

| 代码入口 | 职责/注意事项 |
|---|---|
| `pipeline/step1_generate.py` | PepMLM采样与PPL；历史默认ROOT指向`/home/g/pepmlm_pipeline` |
| `pipeline/step2_bindeval.py` | 早期BindEvaluator结合位点读数；不等同于肽结合概率或当前主reward |
| `pipeline/step3_fold.py`、`step3_graphpep.py`、`step4_analyze.py` | 早期预测、GraphPep与分析入口 |
| `pipeline/feedback_control_1ycr.py`、`feedback_panel.py` | 已完成反馈对照的准备、生成和分析 |
| `pipeline/affinity_diagnostic.py`、`foldx_affinity.py` | 已知亲和力开发诊断 |
| `pipeline/backbone_diagnostic.py`等 | 已完成240结构主链诊断；不默认重跑 |
| `output/receptor_context_mutation_probe_v0_20261022/probe_core.py` | H的无模型核心规则；`probe_run.py`会加载模型 |
| `output/receptor_context_pmi_development_comparison_20261022/paired_core.py` | 实验身份、点效应和误差敏感性等纯CPU规则 |
| `output/receptor_context_disagreement_audit_20261022/audit_core.py` | 配对匹配与分歧分类 |
| `output/bclw_msa_input_audit_20261023/msa_core.py` | 核心投影、多样性与A3M构造 |

部分历史脚本会同步到另一份工作区、下载文件、调用GPU或写既有输出。首次阅读以代码和记录为主，不把`scripts/run_*.sh`或`sync_pipeline.sh`当成当前一键启动命令。

`docs/superpowers/plans/`保存更早的研究实现计划，其中云端或扩展计算提案不代表当前授权或已执行状态；以对应执行记录和本README中的最新状态为准。

## 7. 原工作区环境与路径

| 项目 | 原工作区配置 |
|---|---|
| Windows目录 | `E:\Users\graphpep与pepmlm` |
| 同一目录的WSL挂载 | `/mnt/e/Users/graphpep与pepmlm` |
| 早期另一副本 | `/home/g/pepmlm_pipeline`，与本仓库不是同一份目录 |
| 硬件 | RTX5060 Laptop，约8151 MiB显存 |
| PepMLM/GraphPep Python | `/home/g/miniconda3/envs/graphpep_env/bin/python` |
| ColabFold Python | `/home/g/miniconda3/envs/colabfold_exp0/bin/python` |
| PepMLM权重 | `/home/g/pepmlm_pipeline/models/PepMLM-650M` |
| GraphPep安装 | `/home/g/GraphPep_v1.3` |
| AFM参数 | `/home/g/.cache/colabfold/params/params_model_1_multimer_v3.npz` |
| FoldX授权程序 | 原工作区`tools/foldx/foldx5_Linux/foldx_20261231`，不随Git分发 |

已记录的环境版本包括Python3.11.14、Biopython1.88；ColabFold1.6.3、JAX/JAXlib0.10.2、OpenMM8.6.1。详细证据见[运行环境核对](output/bclw_msa_input_audit_20261023/baseline_runtime_comparability.json)。

`pipeline/colabfold_local_constraints.txt`保存的是更早的JAX0.6.2约束，不代表最新已执行环境；不要根据该旧文件覆盖当前环境。仓库没有提供一键重建完整模型环境的安装包，也不包含权重。换机器运行研究计算前，需要另行恢复资产、路径和许可，并重新核对冻结配置。

## 8. 克隆后可执行的零模型检查

下面四组不加载模型，不预测结构，不下载数据。需要Python3.11或更高版本，在仓库根目录执行。前三组纯规则测试可直接使用Windows Python：

```powershell
python -m unittest discover -s output/receptor_context_pmi_development_comparison_20261022 -p test_paired_core.py -v
python -m unittest discover -s output/receptor_context_disagreement_audit_20261022 -p test_audit.py -v
python -m unittest discover -s output/bclw_msa_input_audit_20261023 -p test_msa_core.py -v
```

H技术原型另外包含两项Linux进程看门狗测试，因此整组需要在WSL执行。原工作区可从PowerShell运行：

```powershell
wsl.exe -d Ubuntu -- /home/g/miniconda3/envs/graphpep_env/bin/python -m unittest discover -s output/receptor_context_mutation_probe_v0_20261022 -p test_probe.py -v
```

在其他机器上则进入WSL中的克隆目录，使用当地可用的Python3.11+运行同一`unittest`命令。看门狗测试只启动短时CPU测试子进程。

本次实际核验：Windows三组共25项通过，WSL原型组8项通过。该结果覆盖上面的33项检查，不代表全部历史pipeline或新机器上的模型推理已验证。

核验随仓库提交的数据字节，在PowerShell运行：

```powershell
@'
import hashlib, json
from pathlib import Path
root = Path.cwd()
manifest = json.loads((root / "docs/repository_data_manifest.json").read_text(encoding="utf-8"))
bad = []
for item in manifest["files"]:
    path = root / item["path"]
    if not path.is_file():
        bad.append(item["path"])
    elif path.stat().st_size != item["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
        bad.append(item["path"])
print({"checked": len(manifest["files"]), "mismatches": bad})
raise SystemExit(bool(bad))
'@ | python -
```

`.gitattributes`使用`* -text`保存原始字节；本地仓库另设置`core.autocrlf=false`。不要对历史冻结文件作批量换行格式化后，再改写原清单使其通过。

旧`audit_msa.py`、`verify_delivery.py`等是当时本地审计入口，可能使用绝对路径、要求335份旧文件（其中含未纳入Git的字节码缓存），或以独占写入保留首次产物。克隆后不要直接重跑它们来覆盖本次结果。仓库数据检查和上述纯规则测试是当前推荐的只读入口。

## 9. 恢复研究时的步骤

1. 先读本README、最新MSA报告和已准备协议；旧交接保留历史，不以其中“未启动”覆盖后来的完成记录。
2. 使用新输出目录，核对本地资产、版本、源代码和权重哈希；原始失败记录与冻结文件保留。
3. 若明确启动已准备的Bcl-w方案，仅按该协议的两seed、预算、参数和停止条件执行；不直接运行硬编码旧kernel目录的runner或summarizer。
4. 结构工程通过之后再单独决定下一项研究。它不会自动补齐Kd、验证H或确定主reward。
5. 亲和力基准仍缺完整数据：PDZ缺合格增强、Homer1融合输入缺项、BH3数量/背景不足，Dutta附件已经读完且未补齐所需Kd/区间。避免重复下载和无上限搜索。

## 10. 本地Git与手动推送

本次建立的是`main`分支的本地研究快照，提交代码、项目说明和数据目录中明确列出的阶段文件。没有配置远程，没有执行push。

推送前可查看：

```powershell
git status --short
git log -1 --oneline
git remote -v
```

由你在GitHub、Gitee等平台创建远程仓库后，在PowerShell执行：

```powershell
$remoteUrl = Read-Host '输入远程仓库URL'
git remote add origin $remoteUrl
git push -u origin main
```

这段命令假设尚无origin且远程可接受当前初始历史；若远程已经有提交，先核对分支关系，不用强推覆盖。后续新增试验目录默认不被Git收录，需要先审查数据范围、在`.gitignore`中明确放行，再更新数据清单。

仓库记录包含原机器绝对路径、运行环境和公开文献数据来源，供溯源使用。第三方软件、模型、结构/数据库记录和文献材料分别遵循各自许可；本次没有为所有内容另行指定统一开源许可证。FoldX程序与模型权重不在提交范围内。

## 11. 方法与数据来源

- PepMLM：[Target sequence-conditioned design of peptide binders using masked language modeling](https://doi.org/10.1038/s41587-025-02761-2)。
- ColabFold：[Making protein folding accessible to all](https://doi.org/10.1038/s41592-022-01488-1)。
- 已知亲和力开发来源：Li等的PMI/p53突变研究，逐项见[实验配对清单](output/receptor_context_disagreement_audit_20261022/experimental_pair_inventory.csv)。
- 公开结构与序列来源：[RCSB PDB](https://www.rcsb.org/)、[UniProt](https://www.uniprot.org/)。
- 突变来源筛查：[SKEMPI 2.0](https://life.bsc.es/pid/skempi2/)，数据库原始记录不自动等于本项目合格标签。
- Bcl-w assay、相关BH3结构与构建体来源，见[BH3来源审计](output/bh3_bclw_source_audit_20261020/report_zh.md)及[结构输入对应表](output/bclw_msa_input_audit_20261023/source_input_correspondence.csv)。

所有“通过”均应结合其验收层次阅读：代码与输入检查通过、结构工程通过、方向评价支持和真实实验改善是不同结论。
