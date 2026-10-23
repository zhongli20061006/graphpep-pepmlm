# Bcl-w 六复合物工程 pilot：配置冻结与运行前核验

日期：2026-10-20（Asia/Shanghai）  
状态：**配置已冻结；未启动预测或评分。** 冻结配置 SHA-256：`6057fa6d5ef225cd979e9756b7760798d745f9e667aa7365670c2db9c20214d0`。

## 当前结论

- 对 Dutta 2013 补充材料的追查未找到图 4c 的逐变体直接 Kd 数值表、误差/置信区间或逐重复曲线；也没有补上满足方向标签要求的增强突变。因此 Dutta 2013 保留为第二个 Bcl-w 突变背景的定性来源，**不计入定量 benchmark 配对数**。补充表中的解离半衰期不能替代 Kd。
- Foight Bcl-w/Bim_A2eT 的 3 条序列、2 个预先固定 seed 共 6 个结构单元，已冻结为**工程可重复性 pilot**。由于候选标签在方法选择中已暴露，这不是盲测，也不验证 scorer 对实验 Kd 的预测能力。
- 六结构最多检查本地输入、预测/放松流程、链映射、预设结构门槛以及 GraphPep/FoldX 等代理读数的跨 seed 稳定性。代理评分与实验亲和力仍分开报告。

## 已完成并核验

- `run_config_frozen_v1.json` 标记 `FROZEN_PRE_SCORE_NOT_RUN`，`run_authorization.status=NOT_STARTED`；没有把配置冻结视为执行许可。
- 输入 FASTA 含 3 个受体:肽复合物条目：受体核心 164 aa、肽 23 aa；序列与冻结来源映射一致。固定比较为每个 seed 内 parent 对两个变体，不跨 seed 混合方向。
- AF2-Multimer 设为本地单序列模式、两个固定 seed、每个输入每 seed 一个模型、3 次 recycle；不请求远端 MSA/模板服务。CPU restrained relaxation 与 GraphPep/FoldX/ESM2 配置、命令、版本、代码/权重/输入哈希均记录在冻结 JSON。
- 结构门槛冻结为：ipTM ≥0.70；肽链 Cα 平均 pLDDT ≥70；双向界面 PAE 最大值 ≤10 Å；链间最小重原子距离 ≥1.5 Å；重原子距离 ≤6 Å 的受体–肽残基对数必须大于 0。身份/质量/接触门槛失败、缺结构、缺分数或非有限分数均保留为失败单元，不替补、不用其他评分器填值。
- 预算固定为最多 6 个复合物单元、累计 GPU 时间 2 小时、该新输出目录新增数据不超过 2 GiB；串行运行。
- 对 10 个冻结输入/实现文件及来源映射、方向判据文件的哈希进行了复核，全部匹配；WSL 指定 GraphPep Python 环境内 ESM wrapper 测试 5/5 通过。该测试不加载真实 ESM checkpoint、不产生实验评分。
- 检查时两个 seed 的结构目录、放松目录和 `scores` 目录都不存在；未发现本 pilot 的预测/评分输出。

## 仍未完成 / 启动前必须检查

1. 启动前即时复查 GPU 可用显存和磁盘余量；当前记录不把旧 GPU 快照当作资源保证。
2. 启动前重新校验冻结配置所列代码、权重、输入与六个确切输出路径；运行中逐项监控 2 小时及 2 GiB 硬上限。
3. 若用户后续明确启动，只运行这个六结构工程 pilot，并保留所有失败项及原始结构；不扩大到全池，不补跑 P27A 最小化。
4. 完整独立亲和力方向 benchmark 仍为 HOLD：当前不够至少 10 个定量合格配对、双向类别和至少两个亲本背景。Dutta 补件未改变这一判定。

## 结果解释冻结

GraphPep 原始分数越低越有利，`parent − variant > 0` 仅指 GraphPep 代理改善；FoldX AnalyseComplex 的实际 A/B Interaction Energy 越低越有利，父代减子代为正也仅指该代理改善。ESM2 masked PLL 是序列先验，接触数与结构质量是描述性结构读数。任何一项都不能写成 Kd 或完整 ΔG(binding) 改善；双 seed 同向也只表示本 pilot 的代理方向稳定。

## 关联文件

- 配置：`run_config_frozen_v1.json`
- 自动核验：`runtime_preflight_v1.json`
- 输入：`complexes.fasta`、`peptide_sequences.fasta`、`graphpep_inputs.csv`
- Dutta 2013 补件审计：`../bh3_second_background_audit_20261020/supplement_followup_zh.md`
- 配置文件清单与 SHA-256：`freeze_manifest_v1.json`、`freeze_manifest_v1.sha256`
