# 独立亲和力评价目标：资料与输入审计工作单

日期：2026-10-07（北京时间）。本工作单采用 writing-plans 的明确文件/验收要求；本阶段为研究准备，直接在本会话完成，不分派代理、不自动进入预测。

目标：交付可评审协议、来源表、原始证据快照、候选与配对清单，并据输入证据判定是否值得进入小规模新预测。

- [x] 完整读取交接及四份最新报告，核对 `E:\AGENTS.md` 和工作区嵌套约束；不以旧 README 覆盖已完成结果。
- [x] 检索原始论文与 PDB，区分 Kd、Ki、IC50、定性 SPOT 与删失数据。
- [x] 在新评分前写出 `frozen_criteria_v1.json`，固定预算、seed、质量/失败口径和逐模型配对汇总。
- [x] 从存档 HTML 提取 Table 2，生成 `candidate_affinity.csv`、`mutation_pairs.csv` 与未筛选的 `literature_table2.csv`。
- [x] 用已安装 Biopython 核对 CIF、UniProt 范围、未解析坐标及 GraphPep 本地训练/测试清单，生成 `input_audit.json`。
- [x] 将旧十肽留在 `development_set.csv`，对主候选核查旧项目精确序列曝光，记录其共享 PMI 骨架。
- [x] 交付 `protocol_draft_zh.md`、`candidate_source_audit.csv`、`decision.json` 及带哈希的准备快照；保留未解决字段。
- [x] 运行只读交付核验，复查来源数值、单位、配对的单点关系、独立性口径、预算与哈希；区分完整性通过和科学进入门槛未通过。

预测、MSA 检索、Amber、FoldX 与 GraphPep 新评分均不是该工作单的默认执行项。通过资料门槛后仍须先冻结精确构建体与运行配置；科学门槛没有通过就保留停止状态。
