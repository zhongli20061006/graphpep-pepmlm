# 独立亲和力评价目标：小型验证协议草案与资料审计

日期：2026-10-07（北京时间）。状态：**评分前判据已冻结；输入尚未放行；本轮新预测 0、新评分 0。**

## 1. 本轮结论与适用范围

已完成交接、四份最新报告、上级约束、原始论文和结构数据库审阅，并准备了候选、配对、构建体、开发曝光和训练清单审计。暂不值得启动新预测：候选的末端化学、SPR Table 2 的误差定义尚未确证；MDMX 实际 assay 构建体还与晶体序列存在待解释差异。

优先候选为 MDM2–PMI-L10A 的 2021 年未使用变体；MDMX–同一肽家族仅作为条件性的同源受体复核。二者不是两个无关家族。它们可支持有限的外部变体检查，但**尚不能满足严格的独立骨架、未富集的增强/减弱验证面板**。论文 Table 2 从 94 条候选中先按 FP 增强筛出 22 条，存在明确选择偏倚。草案中的减弱方向是不同增强突变的反向恢复，不能称为额外测得的失功能样本。

本项目已跑通生成、预测、评分、反馈与一轮重生成；严格 10v10 的 GraphPep 改善各为 2/10，未证明反馈引导优于匹配随机。旧亲和力十肽属于开发面板；FoldX 的总体排序不能消除 P27A 跨 seed 方向不稳定。240 个主链释放结构诊断只有 35/240 通过严格力 RMS 验收，联合 continuation_pass=false 保持不变。本协议不把任何这些结果改称完整自由能验证。

## 2. 原始证据与候选

主要来源为 [Li 等，2021，原始论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC8463443/)（DOI: 10.1016/j.apsb.2021.06.010），Methods 2.1–2.2、Table 2 和 Section 3.1。保留其 HTML 原始快照与完整 Table 2 数值提取，选择过程可回查。

同一 assay 条件：25 °C，SPR 溶液竞争；合成 MDM2(25–109) 或 MDMX(24–108)，分别使用 50 或 100 nM 受体；10 mM HEPES、150 mM NaCl、0.005% P20、pH 7.4；三个独立实验，每个实验技术重复两次。它不是本项目新增的湿实验。

下表的 ± 只表示原文所报误差；当前不把它标成 SD、SEM 或置信区间。每条肽均为 12 aa，表中省略尚未确证的末端化学。

| 肽 / L10A 背景中的变化 | 序列 | MDM2 SPR Kd，nM | MDMX SPR Kd，nM | 草案用途 |
|---|---|---:|---:|---|
| PMI-L10A 基线 | TSFAEYWNLASP | 897 ± 33 | 387 ± 22 | 两受体共用基线 |
| T1L | LSFAEYWNLASP | 123 ± 10 | 91 ± 5 | 两受体减弱反向恢复的起点 |
| A4L | TSFLEYWNLASP | 98 ± 5 | 53 ± 4 | 两受体增强对照 |
| N8A | TSFAEYWALASP | 229 ± 11 | 52 ± 4 | 两受体增强对照；有别于旧 PMI-N8A |
| P12Q | TSFAEYWNLASQ | 132 ± 2 | 134 ± 7 | 仅 MDM2 主清单 |
| S11M | TSFAEYWNLAMP | 198 ± 4 | 106 ± 3 | 仅 MDMX 主清单 |

机器清单只选每受体 5 条，共 10 个受体–肽条目、6 条不同序列。上述完整数值表同时展示未选用受体列，避免隐去来源事实；未选用列不加入主指标分母。

主配对固定如下，每个无向突变边只计一次：

| 体系 | 两个增强方向 | 两个减弱方向（反向恢复） |
|---|---|---|
| MDM2 | L10A → A4L；L10A → N8A | T1L → L10A；P12Q → L10A |
| MDMX，条件性 | L10A → A4L；L10A → N8A | T1L → L10A；S11M → L10A |

每条边都是相同受体、相同 assay、相同 L10A 背景中的单点变化。方向标签用 log10(Kd_parent/Kd_variant)；正为增强，负为减弱。至少三倍变化，并检查 mean ± 2×reported_error 的保守端点是否仍保留三倍方向。这只是误差敏感性检查，**不称为 95% CI**；误差类型和共享基线协方差未明确。八条边的该检查均通过，但这不能替代输入化学核查或证明评分有效。

## 3. 构建体、结构和端化学审计

通过原始 [5UMM](https://www.rcsb.org/structure/5UMM)、[5UML](https://www.rcsb.org/structure/5UML)、3EQY 的 mmCIF 和官方 UniProt 序列逐残基比较：

- 5UMM 的 MDM2 entity 为 85 aa，与 Q00987 的 25–109 范围一致。完整 entity 是 assay 输入的序列核对依据，不能以缺失坐标的 82 aa 链取代它。
- 5UML 的 MDMX entity 为 85 aa，但与 O15151 的 24–108 野生型范围相比有 **Q68A/Q69A/E70A**。3EQY 也有同样三处差异。主文的构建体范围不足以确定 SPR 实际使用 WT 还是 AAA 版本；两种序列均列为待定，不做默认选择。
- 5UMM/5UML 的肽均为 PMI-M3（LTFLEYWAQLMQ），不是本次五条候选。5UMM 的 B 链完整，D 链缺 Q12；5UML 的肽链均缺 Q12 坐标。缺少 OXT 或末端坐标不能证明酰胺化；不能把候选 Kd 贴到 M3 晶体上。
- 论文说明 Boc 合成、HF 裂解与质量验证；访问到的主文未逐条明确本次 SPR 候选 N/C 端化学。自由胺/自由酸只是合成沿革假设，不作为已核实事实。也需记录合成受体的端化学及实际标签。
- Table 3 的 M3 低 pM Kd 来自 ITC，作者还披露多个 ITC 实验中展示低 SD 的代表性一次。本次不将其与 Table 2 SPR Kd、FP Ki 或旧 2010 面板拼接。

原始 CIF、全序列、解析/缺失位置和 OXT 记录在 `structure_source_audit.json`。`receptor_constructs_pending.csv` 和 `candidate_peptides_pending.fasta` 明确标为待定；未生成可直接提交预测的复合物 FASTA。

## 4. 开发曝光与“独立”的边界

旧 2010 PMI/p53 十肽从原始参考 CSV 原样保留到 `development_set.csv`，全部标记 DEVELOPMENT、不能加入验证分母。它们已用于选法和诊断，即使更换 seed 也不会变成盲测。

本轮对工作区 `output` 下所有可识别序列列的历史 CSV 做了精确匹配，六条新序列均无命中；文件清单与哈希见 `input_audit.json`。`candidate_affinity.csv` 中的 `used_for_project_development=False` **仅指该已审计范围未见精确序列曝光**。它不证明所有 WSL 历史、通信、训练语料和同源肽都未曝光。PMI 骨架已被开发，相关 M3 晶体也不是全新家族。

已安装 GraphPep 的训练清单中，`5umm_B` 出现在 `Training_set/list` 第 1837 行；`split_data/train.csv` 从第 5750 行起有 21 行对应记录。Tiam1 备用来源的 `3kze_D` 也有训练命中。本轮仅审计精确 PDB ID 和可访问清单；无命中不等于无同源泄漏，更未证明 AFM、ESM 或 PepMLM 训练无重叠。故不报告 GraphPep 的训练无重叠盲测。FoldX 的有限未使用变体评价另列，不受此标签混淆。

其他原始来源的筛选见 `candidate_source_audit.csv`：pDI/pDIQ 面板主要是 ELISA IC50，不能改写为突变 Kd；Mcl1–Bim 改善主要是 FP Ki 及删失界，构建体/肽端化学不宜跨论文合并；Tiam1 表含预测值、受体突变和多为减弱的肽变化；CAL 化学修饰面板还涉及非标准残基与 FP 表头/方法口径差异。它们均不被填补进本次 Kd 主面板。

## 5. 看新分数前冻结的评价方案

机器判据为 `frozen_criteria_v1.json`，准备快照在 `prepared_snapshot.json`。冻结发生于新预测和新评分均为零时。可补齐原始资料缺项；凡改变候选、端化学假设、主指标、阈值或评分设置，必须先另存版本、注明理由、重新冻结。见分数后不得调阈值、换 seed 或挑选模型救方向；该套数据若用于调参，应转为开发集。

主评价器：现有授权 **FoldX RepairPDB → AnalyseComplex 的 Interaction Energy**。越低越有利，E_parent−E_variant 为代理变化。它不是校准后的 ΔG(binding)、ΔΔG 或 Kd。GraphPep 可作接触/结构的描述性端点，不能用于本面板候选选择或覆盖 FoldX 方向失败；GraphPep 同样越低越有利。

固定 FoldX 配置：pH 7.4、298 K、ionStrength 0.15 M、water=-IGNORE、repair_Interface=ALL、bFactorLimit=0、screen=false。N/C 端参数和实际原子只在化学核实后设置；自由端确认后才能沿用 noNterm=false、noCterm=false。若真实修饰无法被当前流程忠实表达，应停止该面板，不能静默删除修饰。

折叠 seed 为 **20261021 和 20261022**；每 seed 固定 AFM multimer_v3 参数模型 1、2、3，recycles=3、recycle_early_stop_tolerance=0、max_seq=256、max_extra_seq=512、templates=false、unified_memory=false。每次调用只用一个 seed。两 seed 与三模型是计算不确定性样本，不是生物学重复。

每受体使用一份共同受体 MSA，记录来源和哈希并跨变体/seed 复用；肽 MSA 固定单序列。不能让突变各自触发不同肽同源搜索。若使用公共 MSA 服务须记录，不能称为完全离线。输入端化学与运行参数冻结前，不查询 MSA、不加载推理模型。本轮未实现或启动新预测包装器；旧脚本有写死目录，禁止直接执行它们。

统一 CPU Amber 有约束最小化沿用 max_iterations=200、max_outer_iterations=3、tolerance=2.39、stiffness=10、use_gpu=false。保留 AFM 原始结构/置信度、统一 relaxed 结构与 FoldX repaired 结构，不对某个突变选择性增加松弛或释放主链。

结构质量门槛同时要求：ipTM ≥ 0.70；全肽平均 pLDDT ≥ 70；interface PAE median ≤ 10 Å；链间最小重原子距离 ≥ 1.5 Å；至少一对受体–肽残基有 ≤ 6 Å 重原子接触。接口 PAE 的精确定义是：按统一 relaxed 全肽结构找出与肽距离 ≤ 6 Å 的受体残基，再对这些受体残基与**全部肽残基**的两个方向原始 AFM PAE 块取中位数。缺置信度视为未通过；Amber/RepairPDB 不会刷新 AFM 置信度。

## 6. 配对汇总、失败处理和主指标

按 system、parent、variant、同 folding seed、同参数模型配对。每 seed 至少 2 个父/子同时通过质量门槛的模型；先算各配对的 E_parent−E_variant，再取中位数。不能用父与子的独立中位数相减，因为它们可能来自不同模型。

数值死区固定为 0.05 kcal/mol，仅是操作容差，不表示有意义的亲和力效应。每 seed 至少 2 个逐模型差值在预期方向超过死区，且配对中位数也在同方向超过死区，才算该 seed 方向正确。**两个 seed 都正确**才算该边稳定正确。平局、缺配对或方向不一致均为不确定/不正确，保留在分母。

主指标为每受体分别的 `stable_direction_accuracy = 双 seed 稳定正确边数 / 4`。达到 4/4 才能给该有限面板 limited_validation_pass=true；MDMX 若未进入，不记零分，报告“未评估”。两受体分开汇报。共享基线、反向恢复及同源受体存在依赖，不给独立二项置信区间，也不把它们当 8 个独立证据。每家族、每 seed 的 Spearman(proxy, log10 Kd) 仅为次要描述，报 n 和覆盖率；pooled 相关性不能挽救任何方向失败。

无接触、低置信度、明显碰撞、未完整建模、松弛/评分失败均保留原始文件与原因码；没有能量时留空，不能填零或人为有利值。弱结合肽预测失败不自动视为正确减弱。即使工具为无接触结构输出数值，也只放未过滤附录，不进入主方向判断。报告全部 30/60 个预设位置的覆盖情况和所有生成结果，避免只展示通过筛选的模型。

只允许对已识别的技术故障最多一次相同科学参数、相同 seed 重试，保留首轮日志并计入计算时间预算；运行途中未完成不允许补不同 seed。科学方向失败立即结束该面板的验证用途，不能再以新增最小化、P27A 重跑、筛模型或调 FoldX 参数救结果。

## 7. 预算、进入和停止条件

| 项目 | MDM2 主体系 | 条件性 MDMX | 全部上限 |
|---|---:|---:|---:|
| 受体–肽条目 | 5 | 5 | 10 |
| seed × 参数模型 | 2 × 3 | 2 × 3 | 每条 6 |
| 最终预设结构位置 | 30 | 30 | 60 |
| GPU 运行预算 | 4 小时 | 4 小时 | 8 小时 |
| CPU 预算 | 共用 | 共用 | 4 小时 |
| 新目录磁盘 | 共用 | 共用 | 10 GiB |

这是硬预算，不是已测 ETA。每体系先做基线 2 seed × 3 模型的 6 结构 pilot，已包含在 30 内；pilot 墙钟上限 90 分钟。任一 seed 少于两个通过结构门槛的基线模型、资源超限、化学/链序列不一致或未解决的完整性故障，都停止该体系。技术重试不得突破时间/磁盘硬预算；本轮没有 pilot，实际消耗的预测预算为零。

进入顺序与决策：

1. **资料门槛**：确认实际受体全序列/标签、肽和受体端化学、SPR 误差含义与重复层级、源数据来源和开发曝光范围；保持 SPR 内比较，明确增益富集和反向恢复限制。当前未通过。
2. **输入放行**：核实工具可表达真实化学，冻结最终 candidate/pair/complex FASTA、受体 MSA、已有权重和运行版本哈希；新包装器只写新目录，静态检查不引用旧输出。当前仅提供待定序列与工作计划。
3. **预测 pilot**：只在前两门槛通过后启动 MDM2 基线 6 结构。通过才完成余下 24；MDMX 必须先解决 WT/AAA 冲突，再以自身 pilot 决定是否进入。不会默认跑所有文献条目。
4. **有限方向评价**：任何主配对错误或不确定，limited_validation_pass=false；全部 4/4 也只支持本受体、本 assay、所选 PMI-L10A 变体的有限结论。
5. **设计继续门槛**：有限面板通过仍不足以恢复通用闭环优化；需要另一个未富集、含独立测得增强/减弱、最好新骨架的同 assay 数据集，或明确收窄研究主张。旧 P27A continuation_pass=false 不因新面板自动改写。

预算不包括 PepMLM 新生成、全候选池预测、MD/free-energy、重复晶体评分或 P27A 最小化；这些均为零。没有环境重装、GraphPep/FoldX 重下载、付费服务、联系第三方或湿实验订单。

## 8. 本轮已完成、缺什么及下一步

**已完成**：原始来源快照；完整 23 行 Table 2 的单位/assay 保留；10 条候选与 8 条单点配对；旧十肽开发标签；历史 CSV 精确曝光审计；5 个 CIF 与 UniProt 构建体比较；本地 GraphPep 清单重叠审计；判据冻结、未启动任务表及来源/文件哈希；交付完整性核验另见 `verification.json`。

**尚缺、因此不启动**：

- 2021 年所选 SPR 肽实际 N/C 端形式及合成受体端化学的原始支持（例如论文补充资料的逐条合成/质量信息）；不能用一般 Boc 经验补全。
- Table 2 SPR ± 的 SD/SEM/拟合误差定义，以及三个独立实验/每次两重复的汇总口径。
- MDMX assay 用 WT 还是 Q68A/Q69A/E70A 的明确原始记录；晶体范围与 assay 范围相同仍不够。
- 更严格独立验证所需的未富集新骨架同 assay 增强/减弱 Kd 面板；当前反向恢复和共享 PMI 不能充当该证据。

公开补充资料的几次抓取没有得到可用 PDF/压缩包，失败记录明确保留。未把 HTML 验证页、错误响应或未读补充文件当成已核对证据。本轮仅作公开资料审计，不联系作者。后续先补证据，再更新预评分输入版本；有充分支持后，最小的新计算动作才是 MDM2 六结构 pilot。

## 9. 文件索引

- `candidate_source_audit.csv`：两项优先候选与不进入主面板的来源审计。
- `candidate_affinity.csv`、`mutation_pairs.csv`、`literature_table2.csv`：实验标签、配对和完整来源表；不含新代理分数。
- `development_set.csv`：已用于开发的旧十肽。
- `receptor_constructs_pending.csv`、`candidate_peptides_pending.fasta`：待定输入；不可直接启动预测。
- `frozen_criteria_v1.json`、`prepared_snapshot.json`：评分前判据与准备快照。
- `input_audit.json`、`structure_source_audit.json`、`graphpep_overlap_audit.json`：可追溯审计证据。
- `planned_jobs.csv`：60 个未运行预设结构位置，包含可选 MDMX；不是运行结果。
- `source_manifest.csv`、`retrieval_issues.json`：原始来源、派生文本、哈希和失败抓取记录。
- `decision.json`：输入未放行和不启动预测的机器决策。
- `verification.json`、`delivery_manifest.json`：交付完整性与最终文件哈希；不能当作评分有效或生物学验证。
