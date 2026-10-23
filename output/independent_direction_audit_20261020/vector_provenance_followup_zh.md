# Homer1 SUMO 构建体：有限公开追溯

2026-10-20。目的仅为判断已有公开来源能否补齐输入；未登录账户、购买载体、联系作者或运行模型。

Singer 等 2024（[10.1002/pro.5094](https://doi.org/10.1002/pro.5094)）方法4.3报告 pDW363 中表达、体内生物素化的 SUMO–peptide；附件S4提供受体表达序列，但该文和五个附件没有提供相应 SUMO 融合物的完整序列。

追溯同组 Hwang 等 2022（[eLife 70680](https://elifesciences.org/articles/70680)，[PMC8789275](https://pmc.ncbi.nlm.nih.gov/articles/PMC8789275/)）得到补充文件6。该文在 Homer1 论文中被引用，但方法4.3未明确声明所有融合序列、接头和载体克隆沿用2022构建体。2022补充文件6明确给出其 BLI/ITC 大多数肽的表达背景为：

`MAGGLNDIFEAQKIEWHEDTGGSSHHHHHHGSGSGSDSEVNQEAKPEVKPEVKPETHINLKVSDGSSEIFFKIKKTTPLRRLMEAFAKRQGKEMDSLTFLYDGIEIQADQTPEDLDMEDNDIIEAHREQIGGSGSG[PEPTIDE SEQUENCE]`

这是**2022年源文档中的完整前缀**；其中有生物素受体肽、His标签、SUMO以及SGSG接头。可以形成2024样本的有来源重建假设，不能当作已证实的2024克隆序列。尤其不能以同为pDW363或同一研究组，直接消除克隆接头、biotin位点、加工和成熟末端的不确定性。

Addgene 的[pDW363 #8842原始载体记录](https://www.addgene.org/8842/)描述的是pMal-C2来源、用于N端生物素化/MBP融合的基础载体，不能唯一决定后续SUMO改造；该页序列下载需要登录。本轮未尝试绕过访问限制。没有将其他年份或其他蛋白的同名载体代替本次assay样本。

本地证据：`sources/hwang_2022_fulltext.xml`、`sources/hwang_2022_supplementary_files.zip`、`sources/elife-70680-supp6.docx`和`hwang_2022_supp6_text.txt`；下载URL及SHA-256另见`sources/hwang_2022_retrieval.json`和最终产物清单。DOCX按ZIP/XML读取原始文字，本轮未编辑文档。

结论仍为 **HOLD_INPUT**。缺少的是2024实际克隆与该前缀相同的可追溯证据，以及生物素/融合样本到拟评分模型输入的表示约定。它不是“完全没有公开SUMO序列”，也不能写成“已恢复2024精确输入”。若后续选择裁去SUMO/标签、仅预测插入片段，应在任何新评分前另立、冻结并论证截断输入假设，明确评价范围是融合assay标签的方向迁移，不自动等于游离肽亲和力验证。
