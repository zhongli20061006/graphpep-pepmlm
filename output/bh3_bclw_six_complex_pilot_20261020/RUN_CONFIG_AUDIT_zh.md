# Bcl-w 六复合物工程小试：运行配置审计

状态：**输入序列已准备；运行配置仍是草案；没有启动预测或评分。** 这不是亲和力盲测。Foight 2014 候选数据在结构方法和来源审计中属于开发暴露数据；本试点只检查六个 AF2-Multimer 复合物能否按固定输入生成、链映射和双 seed 结构/代理分数是否稳定。

## 已核实

- 工作区中新建的 `complexes.fasta` 有三个条目，每条输入为 Bcl-w 1–164 核心受体加一个 23 aa 肽链，顺序为 receptor:peptide。SHA-256 为 `25c4581aaaac01c1f3ce8ea48a2598bdce982c61502b9a067cc9a91406bf2280`。
- ColabFold 1.6.3、JAX/JAXlib 0.10.2、OpenMM 8.6.1 可导入；`params_model_1_multimer_v3.npz` 已在本地，哈希记入 JSON。当前显卡报告总显存 8151 MiB、空闲 7288 MiB；这只是审计时快照。
- GraphPep checkpoint、ESM2 650M checkpoint 和 FoldX 可执行文件都已存在，哈希及 Python 包版本记入 JSON。不需要重装或重新下载。
- 推荐本地预测用单序列 MSA 模式，避免远程 MSA 服务调用。它偏离此前应用稿的“共享受体 MSA”，因此作为新配置明确记录在首个评分之前；对该小试只解释为工程可行性测试。
- 原始工作区和 `/home/g/pepmlm_pipeline` 的 GraphPep 相关脚本不是同一版本：`feedback_panel.py` 哈希不同，后者缺少当前工作区的 `shared_receptor_msa.py` 和 `relax_feedback_panel.py`。没有执行同步。运行时应从 `/mnt/e/Users/graphpep与pepmlm/pipeline` 调用当前工作区代码，并单独校验 GraphPep 安装与兼容 mapper。

## 开跑前还需冻结

1. ESM2 的 masked pseudo-log-likelihood 公式已写入配置草案，但对应可执行 wrapper 和 pooling 实现尚未冻结；应先固定代码哈希，再看任何新分数。
2. 原方案只写“2 GiB”，没有指定资源类型。草案建议按新实验输出目录的磁盘上限解释；同时 GPU 时间上限为 2 小时。超过任一预算即停止，不重抽失败结构。
3. 真实预测/评分是新的计算实验，尚未启动。其余准备工作是可审阅的本地文件与只读环境核验。

## 解释边界

源论文的亲本、I3dL 增益候选与 I2dM 减弱候选来自同一 Bcl-w 直接 Kd assay，但该系统仍只有一个工程肽背景，不能补足冻结完整面板的最少样本/背景门槛。GraphPep 的 parent−variant 为正只代表代理分数改善；FoldX Interaction Energy 也不是完整 ΔG(binding)。ESM2 与接触数只作独立描述性 baseline。双 seed 同向也仅表示这个小试的结构/代理稳定，不证明 reward 有效。

机器配置见同目录的 `run_config_draft.json`；冻结输入映射仍以来源审计目录中的 `input_representation_frozen_pilot_v1.json` 为准。