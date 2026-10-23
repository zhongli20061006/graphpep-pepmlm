# 亲本快速内核路径诊断：准备稿已冻结，未启动

协议：`BCLW_PARENT_KERNEL_ABLATION_V1`。可执行参数及输入哈希以 [协议JSON](next_compute_protocol_prepared_v1.json) 为准。它不预指定任何主reward。

**明确问题：** 在同一受体核心/亲本序列、两个既定seed、同一模型与回收设置下，关闭快速内核路径是否足以恢复亲本结构置信度？这不是“关闭就一定修好”的预判，也不验证亲和力。

| 项目 | 冻结内容 |
|---|---|
| 输入 | `prepared_inputs/parent_only.fasta`：Q92843 1–164核心 + Bim_A2eT标准23-mer；仍省略染料、酰胺、c-myc/未知接头；不改Q133 |
| 基线 | 已完成pilot中的两个同seed亲本；不重跑、不重评分 |
| 新任务 | 仅亲本×20261008/20261009，最多2次预测；同一CPU restrained Amber设置，最多2个放松结构 |
| 唯一预测选项改变 | `--use-fast-kernels` → `--no-use-fast-kernels`；输入路径/输出目录因仅亲本及新试验相应变化 |
| 固定项 | 已安装ColabFold1.6.3、Multimer v3 model1、3 recycle、early stop tolerance0、ensemble1、single_sequence、无模板、unpaired、max-seq256/max-extra512（生效1/1）、compile-mode fast |
| 环境 | 现有colabfold_exp0；现有model1参数；启动时核对版本/源码/权重哈希，不升级/重装/下载 |
| 预算 | GPU标记进程累计wall≤1800秒；CPU放松wall≤900秒；输出≤536870912字节；串行；必须有预算监测 |
| 输出 | 新的`output/bclw_parent_kernel_ablation_v1`；若已占用，先冻结唯一新目录，不覆盖旧结果 |
| 明确排除 | 所有突变预测、GraphPep/FoldX/ESM2、PepMLM、远端MSA/模板、MD/主链最小化、训练、付费与湿实验 |

原结构准入阈值保持：链A164/B23精确对应、非零接触、ipTM≥0.70、肽平均Cα pLDDT≥70、全双向跨链PAE最大值≤10 Å、相同放松后最小链间重原子距离≥1.5 Å。额外要求受体均值pLDDT≥70，只作为新诊断准入；原pilot仍是0/6，不回溯改写。

1. **2/2全过：** 只支持这一实现路径改变足以在两个seed下恢复亲本；不能据此宣布reward有效。突变与评分仍需另立具体方案。
2. **1/2全过：** seed敏感，保持HOLD，不再抽seed。
3. **0/2全过：** 该选项改变单独不足；不能推出所有环境无问题或MSA/端基必然是原因。
4. 只改善均值但未全过，仍是失败。OOM、NaN、程序错误或预算触顶时停止，保留产物与失败记录，不自动降级、换模型、加recycle或重试。

未来启动前检查：适用AGENTS及最新指令、原文件完整性、输入及参数哈希、现有GPU进程/空闲显存、磁盘余量、预算监测和输出目录。CLI帮助已经确认支持`--no-use-fast-kernels`；本轮没有运行CLI预测或加载参数来测试此分支，所以RTX5060约8GB上的该分支可运行性尚未验证。

**授权状态：本轮只准备，不启动。** 未来用户明确要求执行本已准备协议，即按该范围启动，无需再次确认相同许可。

另一候选是共同受体MSA分支。当前没有可冻结的Bcl-w同源MSA；现有六个A3M只是query序列，不能拿重复query冒充同源。不得把`single_sequence`标志原样用于有MSA的A3M后声称检验MSA：本机代码会在此模式重置为query-only。待输入来源、查询一致性、覆盖/重复/污染与实际深度审计后，再独立冻结仅改变MSA的方案。给定实验pose的分支也另立协议，明确其科学问题变化。
