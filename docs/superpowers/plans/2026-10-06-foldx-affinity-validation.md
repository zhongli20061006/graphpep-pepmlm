# FoldX MDM2 affinity validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. Track execution with checkbox steps.

**Goal:** 用 FoldX AnalyseComplex 的链间 Interaction Energy 检验已有十条 MDM2 肽的实测 Kd 方向，先冻结参数与判据，不重新生成肽或折叠。

**Architecture:** 原 affinity diagnostic 数据保持不变；新目录 `output/foldx_affinity_mdm2_20261006` 保存 60 份重原子结构、配置、冻结协议、真实 FoldX 原始输出及分析。工具暂未安装：用户正在获取官方授权 Linux 包；先完成可独立核验的准备工作，缺工具时不得输出预测能量或通过结论。

**Tech Stack:** 原 WSL graphpep_env Python、标准库、既有 pandas/SciPy；用户提供的官方 FoldX Linux executable。无新模型训练，无新 GPU 作业，无 Git 元数据。

## Task 1: 输入与协议准备

**Files:** `pipeline/foldx_affinity.py`, `tests/test_foldx_affinity.py`, `output/foldx_affinity_mdm2_20261006/`。

- [x] 核对旧 `reference_affinity.csv`、`protocol.json`、60 行 `analysis/model_results.csv`，逐行检查结构 SHA256；十条序列、两个 seed、三个模型身份必须完整且唯一。
- [x] 保留每个原 AFM PDB 的重原子坐标，移除 H/D，不移动坐标、不改序列。每个输入均为 A 链 85 aa、B 链 12 aa，残基连续，无重复重原子，骨架 N/CA/C/O 完整。
- [x] 每个模型建立独立工作目录，写 `input.pdb`、`repair.cfg`、`analyse.cfg`；模型名使用 ASCII，以避免 FoldX 对 Unicode 路径解析不确定。
- [x] 冻结 RepairPDB 后 AnalyseComplex；分析 A,B；temperature=298 K、pH=7.4、ionStrength=0.15 M、water=-IGNORE、screen=false、noNterm=false、noCterm=false，保持游离末端处理。298 K 是对文献室温的统一近似，0.15 M 对应 NaCl 的主要贡献，并非精确复现全部 assay 溶液。
- [x] 保存源数据、源结构、导出结构、配置及本轮脚本 SHA256；重跑 prepare 只核验冻结内容，不静默覆盖。
- [ ] 工具取得后另存 immutable tool manifest（可执行文件、帮助/版本输出、配套资源哈希）；FoldX 5 内嵌参数在调用前由 executable 哈希绑定，自动生成的 rotabase 文件在实际生成后逐作业记录哈希。

## Task 2: 真实评分与恢复

**Files:** `scripts/run_foldx_affinity.sh`, `tools/foldx/README.md`。

- [x] 脚本使用 `/home/g/miniconda3/envs/graphpep_env/bin/python`；`prepare` 不依赖 FoldX；`run --foldx /absolute/path` 在缺 executable 时清楚报错。
- [ ] 先处理两个父代以及 N8A、P27A、L22A，再处理其余，固定为全部 60 个。每个 FoldX command 保留 stdout/stderr 与命令、耗时、返回值；失败记录为失败，不选择性丢弃弱结合对照。
- [ ] RepairPDB 使用独立配置，在原 source 文件外工作；验证 repaired PDB 序列仍为 85+12，再 AnalyseComplex。只从含 `Interaction Energy` 的原始 interaction/summary 表读取对应 A/B 项；不能读取 BuildModel stability Dif_。
- [ ] 原始输出解析要求唯一匹配行、链对正确、能量有限；科学数据不得从终端欢迎信息或文件首行任意数字提取。恢复缓存以输入/配置/脚本/tool 指纹与 raw/repaired 文件哈希为依据。

## Task 3: 固定分析与独立验证

**Files:** 输出 `analysis/model_results.csv`, `sequence_summary.csv`, `mutation_directions.csv`, `rank_correlations.csv`, `report_zh.md`；独立 verifier。

- [ ] 每个 seed 内按模型号配对：gain = median(parent energy − variant energy)，正值改善。主质量过滤每个 seed 至少两模型；全模型统计亦报告。不把结构模型当作生物学重复。
- [ ] 七项主要方向 ≥3 倍实测变化：N8A/P27A 变强，PMI Y6A/F3A/W7A 与 p53 L22A/L26A 变弱。K24A 描述性，不纳入主方向数。
- [ ] 沿用门槛，每 seed 十条肽全部可比较、总体 Spearman(energy,log10Kd)≥0.6、PMI/p53 家族各>0、至少6/7正确且两项增强正确、七项全部可比较。分数波动、无法区分与失败均如实报告，不改变门槛求通过。
- [ ] 只有固定 60 个作业全部完成尝试（成功或明确失败）才允许总体继续判据为 true；40/60 等运行中快照只能作暂时统计，不得提前宣布通过。
- [ ] 对真实结果独立从 raw 重新解析与复算，核对 SHA256、chain、单位、paired median 与恢复一致性。科学 guard 测试覆盖 parser 拒绝非有限/歧义/错链、配对统计，另用明确合成数据检查批处理接口，不能把合成输出列为真实 benchmark。
- [ ] 若 FoldX 不通过，不启动当前目标的新子代生成。MD/MMGBSA 仅为后续候选且仍须独立检验；不能用一次最小化或短 MD 声称自由能收敛。
- [ ] 若通过，冻结选定方法，再加入未用于选方法/调参的已知突变面板；训练集重叠未审计时不称严格外部测试。新子代实验 Kd 阶段需要实际实验资源，此次计算不代替实验。

## Ready-to-run commands

```bash
cd /mnt/e/Users/graphpep与pepmlm
/home/g/miniconda3/envs/graphpep_env/bin/python pipeline/foldx_affinity.py prepare
/home/g/miniconda3/envs/graphpep_env/bin/python tests/test_foldx_affinity.py
bash scripts/run_foldx_affinity.sh /absolute/path/to/authorized/foldx
```

## External dependency status

FoldX PATH、WSL /home/g、/opt、/usr/local 和常用 Windows 目录暂未找到软件。官网要求注册登录后下载；用户已回复“还没有，我去官网下载”。用户下载期间完成输入准备与测试，实际 FoldX 评分保持 pending。

## Verified preparation results

- 真实 60 个模型输入及 120 个配置已导出；独立输入核验通过，源结构未变、所有重原子坐标保持一致、完整 85+12-aa 序列及游离末端保留；移除原 Amber H/D 共 50,214 个。
- 直接 WSL Python unittest 共 8 项通过，最近一次 17.991 秒，exit 0。合成原始表、失败保留、恢复与统计核验只在 TemporaryDirectory 内测试；不属于真实 FoldX 评分。
- 科学审查与代码复审通过；缺失 K24A 不允许用九条肽宣称通过，40/60 等部分运行不允许提前继续。
- 真实 source/input/script/config SHA256 已冻结，独立输入审计与重复核验通过。真实能量 `ok=0, failed=0, not_run=60`；两个 seed 均 `NOT_EVALUABLE`，不是新指标验证失败。
- Task 2 的实际 RepairPDB/AnalyseComplex 和 Task 3 的真实统计/评分审计、后续 holdout 尚未执行；依赖用户提供官方授权 executable。
