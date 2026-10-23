# Bcl-w 亲本快速内核分支结果

协议：`BCLW_PARENT_KERNEL_ABLATION_V1`。2个固定seed各完成1次预测和1个同设置CPU放松；未调用GraphPep、FoldX、ESM2或其他reward。

| Seed | ipTM | 受体均值pLDDT | 肽均值pLDDT | 双向界面最大PAE (Å) | 原始→放松最小链间距离 (Å) | 放松后接触残基对 | 全门通过 |
|---:|---:|---:|---:|---:|---:|---:|:---:|
| 20261008 | 0.12 | 40.49 | 36.45 | 30.09 | 0.496 → 2.593 | 46 | 失败 |
| 20261009 | 0.25 | 40.92 | 42.28 | 29.86 | 0.778 → 2.585 | 38 | 失败 |

通过数为 0/2。冻结分支结论：HOLD; disabling fast kernels alone did not recover the parent under these two seeds; this does not establish that other runtime or input factors are innocent or causal.

生效配置复核：除输入从同seed的三条件面板缩为亲本单条件（`num_queries: 3→1`）及指定的 `use_fast_kernels: true→false` 外，没有其他配置差异。有效MSA仍为single-sequence（每条链1行）、无模板、Multimer v3 model 1、3 recycle。 两个新生成A3M均核对为 `#164,23	1,1` 的单序列分链填充格式，没有同源序列；四个recycle步骤（0–3）均出现在日志。细节见 `input_artifact_verification.json`。

两组原始结构与放松后结构的A164/B23序列和连续编号均已核对。原始结构存在小于1.5 Å的链间近接，但按冻结判据使用相同CPU放松后几何；若放松后距离通过且接触非零，几何门通过。

此分支只回答快速内核选项是否足以改善这两个固定seed下的结构置信度。结构门失败不等于实验不结合，也不定位具体根因；本次没有新Kd、完整结合ΔG或亲和力方向准确率。

复核了旧pilot manifest中的265项，哈希仍全部一致。实际GPU/CPU子进程耗时、磁盘和原始结果哈希见 `budget_reconciliation.json`、`execution_record.json` 与 `artifact_manifest.json`。

## 文件

- `quality_gate_summary.csv` / `.json`：冻结门槛逐项结果及同seed基线对照。
- `input_artifact_verification.json`：新生成A3M、有效配置、日志回收步骤与产物数复核。
- `decision.json`：按2/2、1/2、0/2规则冻结的结论。
- `budget_reconciliation.json`：按各子进程实际wall time核算的预算记录。
- `execution_record.json`、`logs/`、逐seed结构目录：原始执行与产物。
