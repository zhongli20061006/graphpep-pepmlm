# Bcl-w 六复合物工程 pilot 运行结果

日期：2026-10-20（Asia/Shanghai）  
状态：**六个结构预测及预设代理评分流程已执行；0/6 结构通过冻结质量门槛，因此没有可解释的突变方向结果。** 本运行不是独立亲和力验证，也不是盲测。

## 运行范围与资源

- 输入是 3 条 Foight Bcl-w/Bim_A2eT 肽序列，每个固定 seed 各预测一次，共 6 个 AF2-Multimer 复合物；两个 seed 均按冻结配置使用本地 single-sequence 模式。六个原始结构、置信度 JSON 和 restrained Amber 放松结构均已保存。
- GraphPep 6/6、FoldX RepairPDB→AnalyseComplex 6/6 成功产出原始读数；ESM2 masked PLL 为 3 条序列各算一次。一次 ESM 初始加载因 PyTorch 安全反序列化兼容问题失败，随后用单独冻结的兼容入口重试成功；第一次失败日志保留。兼容入口仅在 `weights_only` 安全加载中允许该检查点实际含有的 `argparse.Namespace`，没有改用 `weights_only=False`。PyTorch 文档将 `safe_globals` 作为 `weights_only` 加载的受限白名单机制。[官方序列化说明](https://docs.pytorch.org/docs/stable/notes/serialization.html)
- 累计 GPU 计量时间 213.3 秒（约 3.6 分钟），低于 2 小时上限；放松与 FoldX CPU 阶段共约 533 秒。pilot 目录新增数据约 11 MB，低于 2 GiB 上限。没有远端 MSA、模板服务、训练、全池预测或 P27A 最小化。

## 结构身份和质量门槛

6/6 结构的链 A/B、164 aa 受体与 23 aa 肽序列及连续残基编号均匹配；6/6 有非零接触，最小链间重原子距离为 2.405–2.606 Å，没有触发 1.5 Å clash 门槛。但所有结构均同时未达到 ipTM、肽链平均 Cα pLDDT 和界面 PAE 门槛：ipTM 为 0.11–0.26（要求 ≥0.70），肽链 pLDDT 为 36.01–42.49（要求 ≥70），界面 PAE 最大值为 29.75–30.28 Å（要求 ≤10 Å）。因此 0/6 进入方向可评估集合。

| Seed | 复合物 | ipTM | 肽链平均 pLDDT | 最大界面 PAE (Å) | 最小链间距离 (Å) | 6 Å 接触残基对 | 质量门槛 |
|---|---|---:|---:|---:|---:|---:|---|
| 20261008 | parent | 0.12 | 36.44 | 29.78 | 2.405 | 53 | 失败 |
| 20261008 | I3dL gain | 0.20 | 39.44 | 29.84 | 2.485 | 39 | 失败 |
| 20261008 | I2dM loss | 0.11 | 36.01 | 29.75 | 2.423 | 63 | 失败 |
| 20261009 | parent | 0.26 | 42.49 | 29.86 | 2.470 | 32 | 失败 |
| 20261009 | I3dL gain | 0.26 | 41.39 | 30.14 | 2.606 | 29 | 失败 |
| 20261009 | I2dM loss | 0.20 | 40.82 | 30.28 | 2.543 | 40 | 失败 |

## 代理分数：仅保留原始诊断读数

配对差值只在同 seed 内计算。GraphPep 越低越有利，表中为 `parent − variant`；FoldX 为实际 A/B Interaction Energy（kcal/mol），也是 `parent − variant`。正值只表示对应代理读数更有利。尽管评分程序成功运行，相关结构未通过质量门槛，所以下表所有配对均标记为**不可评估**，不判为成功、失败的亲和力方向或 seed 稳定改善。

| Seed | 变体 | GraphPep parent−variant | FoldX parent−variant (kcal/mol) | 冻结方向判定 |
|---|---|---:|---:|---|
| 20261008 | I3dL gain | −0.0655 | +2.9286 | 不可评估：结构门槛失败 |
| 20261009 | I3dL gain | −0.1178 | −2.7314 | 不可评估：结构门槛失败 |
| 20261008 | I2dM loss | −0.4672 | +5.4655 | 不可评估：结构门槛失败 |
| 20261009 | I2dM loss | −0.0093 | +1.9386 | 不可评估：结构门槛失败 |

FoldX 对 I3dL 的原始差值跨 seed 变号；对 I2dM 两个 seed 都为正，但该变体在来源标签中是减弱候选。因为模型结构全部失败，不能据此解释为模型反例或真实亲和力趋势。没有计算 pooled 相关性。

ESM2 是序列先验，不是界面或亲和力评分。相对 parent 的总 masked PLL 差值为：I3dL `+0.4447`，I2dM `−0.3958`；正值仅表示该模型对变体序列的先验偏好较高。跨 seed 结构描述读数为：parent 肽主链 Cα RMSD 3.72 Å、接触 Jaccard 0.518；I3dL 分别为 39.70 Å、0.619；I2dM 分别为 5.37 Å、0.392。它们不改变结构质量失败状态，也不表示结合变化。

## 判定与后续

这轮完成了预测、放松、链映射审计和代理执行，但在最前置的结构置信度关口整体失败。**停止把这组分数用于奖励函数比较、方向判断或亲和力结论；不改阈值、不换 seed 重跑。** 下一步应先审查当前单序列 AF2 输入与 Bcl-w/Bim 体系的结构可识别性，再决定是否另立一种结构输入/已知复合物模板方案；这需要独立冻结，不能把本轮失败结构分数当作验证结果。完整独立亲和力 benchmark 仍为 HOLD，Dutta 2013 的补充材料也没有补齐定量双向标签。[Foight 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4520770/)、[Dutta 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3631442/)

## 文件

- 冻结配置与预评分审计仍在 `run_config_frozen_v1.json`、`launch_preflight_v1.json`；预评分冻结记录未改写。
- 结构门槛行级结果：`structure_gate_audit_v1.csv` 和 `.json`。
- 同 seed 代理配对及序列先验：`paired_proxy_summary_v1.json`、`paired_proxy_summary_v1.csv`。
- 原始代理表：`scores/graphpep.csv`、`scores/foldx.csv`、`scores/esm_pll.csv`。
- FoldX 每个复合物的 RepairPDB/AnalyseComplex 配置、原始输出和哈希：`foldx_jobs/`。
- 运行命令、耗时、失败重试和预算：`run_execution_record_v1.json`、`logs/`。
- ESM 安全加载兼容补充：`esm_loader_compat_addendum_v1.json`、`esm_pll_pytorch27.py`。
- 全部结果文件哈希：`run_result_manifest_v1.json` 和 detached `.sha256`。
