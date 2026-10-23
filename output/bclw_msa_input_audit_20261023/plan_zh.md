# Bcl-w MSA输入审计执行计划

日期：2026-10-23，Asia/Shanghai。用户已授权执行此前提出的输入审计与条件方案准备。

**目标：** 完成已有93条UniProt候选的可复算筛选、对齐和A3M输入检查，按原冻结门槛交付GO/HOLD。
**实现：** 新目录内独立CPU脚本，使用已安装Biopython 1.88；仅读取原始记录、已安装解析器源码及历史结果。
**授权边界：** 本轮预测、放松、PepMLM/GraphPep/FoldX/ESM2打分、训练、下载和外部服务调用均为0。输入合格仅准备未来亲本双seed诊断，不能执行。
**预算：** 自本次任务开始约1小时人工/工具审计；CPU工作进程限300秒，新增输出限16 MiB。超限保留失败，不改变筛选门。

## 文件与任务

- [x] 冻结剩余实现约定，记录旧来源与冻结文件的SHA256；原目录不写入。
- [x] `test_msa_core.py`先用合成序列验证插入、缺失、164/23链边界、去重与多样性计算；运行后再实现`msa_core.py`。
- [x] `audit_msa.py`读取原始93条记录，逐条输出排除理由、全长比对、核心投影、最终代表序列和多样性。
- [x] 核对原始Q92843全长、1–164核心、Q133，以及旧亲本FASTA。冻结标准肽与省略化学保持原样。
- [x] 输出受体A3M及亲本复合物A3M，用已安装解析器的函数体作CPU隔离解析核对；不导入ColabFold/JAX模型运行时。
- [x] 按10月22日原冻结门判断输入；任何门失败则HOLD并保留候选输入，不扩来源或放宽阈值。
- [x] 仅在输入门及格式检查通过时写出未来双seed诊断协议，包括输入、命令、必要MSA参数差异、预算及停止条件。
- [x] 核验新产物与旧目录完整性，交付中文报告与机器可读decision。

## 实现约定

原规则完整读取`../bclw_construct_msa_audit_20261022/scope_and_selection_frozen.json`，本计划不替代该文件。
新增约定仅补齐可复算的排序、身份分母、聚类与插入边界，详见`selection_implementation_frozen.json`；在真实序列比对前冻结。

固定执行顺序：原始字段/长度/片段/融合/标准残基检查 → 全长global alignment → 核心覆盖/身份检查 → 排除query同核心 → 每taxon选一条 → 核心列去重 → 最多64条代表选择 → 输入门。

全局比对采用原冻结BLOSUM62、open=-10、extend=-0.5；对Q92843全长193 aa对齐后仅投影1–164，不将全长疏水尾段带入核心。
该MSA是有来源记录的BCL2L2基因注释同源候选比对；不声称完整进化搜索、严格正交关系或结构预测一定恢复。

## 运行与验收

使用WSL当前挂载工作区，不使用`/home/g/pepmlm_pipeline`副本。

```text
/home/g/miniconda3/envs/graphpep_env/bin/python -m unittest discover -s output/bclw_msa_input_audit_20261023 -p test_msa_core.py -v
/home/g/miniconda3/envs/graphpep_env/bin/python output/bclw_msa_input_audit_20261023/audit_msa.py
```

验收需要：93项去向可核对；A3M去小写后每条核心164列、复合物187列；链序列与旧输入一致；无伪造肽同源配对；独立复算指标与门判一致；旧文件哈希未改变。
测试验证输入处理，不构成生物学验证。工程输入GO也不释放正式Kd标签，不证明H或任何reward有效。

工作区不是git仓库；交付留在新输出目录，不新建仓库、不提交旧文件。
