# 已有H/d分歧与本地实验配对：有限只读审计

用户授权“行，继续尝试吧”，承接的是从已有读数提取H与d分歧、查本地实验记录。本轮零模型前向、零训练、零下载、零结构预测、零新增GraphPep/FoldX/ESM2调用；只新写本目录。

执行前固定：沿用1e-4 nats数值容差与PMI开发对照的方向规则。d按正负判断；H要求两固定控制修正同向，否则拒绝。严格反向（gain/loss对立）、H拒绝、d平局后H有方向分别统计，拒绝不能冒充方向反转。raw H符号仅作辅助诊断。

合并首轮60输入probe和随后19输入开发对照的base分布，按checkpoint、完整受体序列、完整亲本肽、位点、FP32/eager/TF32off去重。重复和复用不计新位点。预期19个位点（p53 13mer 8个、PMI 12mer 11个），每位点19个非WT替换，共361行；同时标记原283个生成候选池与后来仅用于评价的位点/历史排除，绝不将后者提升为推荐。

实验核对范围预先限定：

1. Li 2010原文两表的既有source_input_audit.csv，完整亲本与变体双端精确匹配。保留近中性、ND、p53输入/重复数缺项，全部formal标签仍HOLD。
2. Li 2021本地literature_table2.csv内完整22个MDM2/PMI-L10A原始forward配对，candidate_affinity.csv提供已核对的构建体/端基/误差元数据。22个是FP富集后测Kd的集合，不保证每个SPR点效应均达两倍。亲本L10A不能冒充普通PMI，不使用旧mutation_pairs.csv反向边凑减弱。
3. 只查既有SKEMPI索引已确定的3EQS_A_B、3LNZ_A_B、3EQY_A_C对应原始行。前两者属于同一Li 2010数据，其中3LNZ为N8A反向边；3EQY是MDMX。只作来源去重/受体排除记录，不算新独立标签，不重新筛全库找lead。

禁止仅按突变名或variant序列拼接，不能裁13mer来继承12mer标签，不能换受体/assay/端基。近似匹配只能进入单独的HOLD表。

结果按分数找分歧再查已见标签，明确为后验开发诊断，不是独立测试。完整361行与分歧子集都保留，避免把子集成功率当总体性能。即使无合格分歧也如实交付，不新增计算、放宽两倍效应或借点估计猜实验方向。

预计输出：all_predictions.csv、disagreements.csv、experimental_pair_inventory.csv、exact_experimental_matches.csv、near_matches_hold.csv、skempi_named_lookup.csv、source_scope_audit.csv、summary.json、report_zh.md、verification.json、不可变输入/代码/交付清单。CPU进程<=120秒，目录<=16MiB；无GPU预算。
