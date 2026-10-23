# 执行与解释边界

- 原始亲和力数据由独立核对者逐条核实为 Li et al. 2010 Table 1/3 的 MDM2 数值。当前基准十条肽均为化学序列明确的 12-mer；无 acetylation/amidation，不把 1YCR 解析到的 13 个残基当作 p53(17–28) 实验肽。
- 3LNZ 保留八个 receptor–peptide 副本，几何配对与 deposited A–B/C–D/E–F/G–H/I–J/K–L/M–N/O–P 一致。共同已解析范围为受体 MDM2 27–108 与肽 1–11；原始 native 重原子坐标不松弛或人为补全。
- 输入准备期间修复作者编号映射、裁剪结构缓存字段继承及评分依赖溯源。AFM 预测开始后、任何 GraphPep 评分或亲和力分析前，补齐有限值质量检查、独立 ESM 嵌入缓存和实际裁剪接触计数。初始与中间协议均保留。`protocol_revision_consistency.json` 独立核对三版旧协议的十个 MSA、亲和力 CSV、native manifest 哈希以及 seed/判据，全部一致；没有依据结果换样、换种子或改成功门槛。
- 四项针对单位/序列、相关方向、配对差值中位数和跨结构缓存污染的检查通过；使用现有 Python 直接执行测试，没有为缺失 pytest 安装新依赖。
- GraphPep 源码、本地兼容补丁、atom-pair potential、评分 checkpoint 与 ESM checkpoint 均纳入冻结哈希；本轮独立 embedding cache 防止之前的无版本缓存影响基准。
- 裁剪 AFM 记录中的 ipTM/pLDDT/PAE 与 `structure_pass` 来源于完整 85+12-aa 预测，裁剪不会重新计算模型置信度。裁剪评分使用 82+11-aa 的实际坐标、接触计数与序列；实测 Kd 属于原完整 12-mer，不能解释为截短肽的 Kd。
- 两个 seed 和三个 AFM 参数模型用于检查计算敏感性，八个晶体副本用于检查构象敏感性，均不能计作独立亲和力测量。总体相关系数可能被破坏核心锚点的巨大变化主导，须同时检查家族内排名和局部突变方向。
- 本基准未审计 AFM/GraphPep 训练集重叠，不能称为严格的独立外部测试。若此小基准失效，可以指出当前目标在这些已知对照上的限制；若通过，也不能证明跨受体、跨肽骨架泛化。
- 文献 Kd 是现有已知肽的外部标尺；本轮没有测量之前生成子代的 Kd，没有计算 ΔG(binding)，Amber 仅用于一致的几何最小化。
