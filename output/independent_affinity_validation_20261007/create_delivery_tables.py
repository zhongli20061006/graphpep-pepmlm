"""Write research audit tables and a source manifest; no model/score execution."""
from __future__ import annotations

import csv
import hashlib
import html
from pathlib import Path
import re


TASK = Path(__file__).resolve().parent


def write_csv(name: str, rows: list[dict[str, object]]) -> None:
    path = TASK / name
    if path.exists():
        raise FileExistsError(f"Refuse overwrite: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit_row(**values: object) -> dict[str, object]:
    fields = ["system_id", "priority", "source_doi", "source_url", "source_location",
              "receptor_construct", "construct_tags_or_mutations", "peptide_length",
              "peptide_N_terminal", "peptide_C_terminal", "assay_and_conditions",
              "quantity_reported", "affinity_evidence", "error_and_replicates",
              "enhancement_and_loss_evidence", "PDB", "structure_vs_assay_match",
              "used_in_current_project_development", "training_overlap_audit",
              "selection_bias_or_independence", "missing_information", "current_decision"]
    assert set(values) == set(fields), set(values).symmetric_difference(fields)
    return {field: values[field] for field in fields}


common = {
    "source_doi": "10.1016/j.apsb.2021.06.010",
    "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8463443/",
    "source_location": "Methods 2.1/2.2; Table 2; Section 3.1; local sources/PMC8463443.html",
    "peptide_length": "12 aa; 每体系5条; 两受体合计6条不同序列",
    "peptide_N_terminal": "未确证; Boc/HF 合成不作为自由胺已核实证据",
    "peptide_C_terminal": "未确证; 缺失 OXT 不能作为酰胺证据",
    "assay_and_conditions": "SPR solution competition; 25 C; 10mM HEPES/150mM NaCl/0.005%P20/pH7.4",
    "quantity_reported": "Kd nM; 与同表FP Ki及Table3 ITC Kd分开",
    "error_and_replicates": "3 independent experiments each duplicate; ± 类型与汇总口径待核; 不填SD/SEM",
    "used_in_current_project_development": "六条所选精确序列在已审计工作区历史CSV无命中; PMI家族已用于选法; 非完整历史排除证明",
    "selection_bias_or_independence": "94条先经FP增强富集选22条; 共享开发PMI骨架; 反向恢复不等于新增失功能样本; 非盲测",
}
rows = [
    audit_row(**common,
              system_id="MDM2_PMI_L10A_2021", priority="1 主候选但未放行",
              receptor_construct="human synthetic MDM2(25-109); 85 aa",
              construct_tags_or_mutations="5UMM entity1与Q00987 25-109完全一致; 主文未述标签; 端化学/逐试剂身份尚不完整",
              affinity_evidence="SPR Kd nM: L10A897±33; T1L123±10; A4L98±5; N8A229±11; P12Q132±2",
              enhancement_and_loss_evidence="增强 L10A→A4L/N8A; 减弱 T1L/P12Q→L10A 为两个不同增强边反向; 每边仅计一次",
              PDB="5UMM 1.65A; M3 peptide LTFLEYWAQLMQ; 3EQS为相关旧PMI结构",
              structure_vs_assay_match="受体范围/全序列一致; 肽不一致; 所选候选无精确PDB已确认; 不将其Kd赋给M3晶体",
              training_overlap_audit="5umm_B: GraphPep Training_set/list line1837; train.csv firstline5750 21rows; AFM/ESM训练清洁未证明",
              missing_information="所选SPR肽及合成受体端化学; SPR误差定义; 更严格独立新骨架未富集增强/减弱面板",
              current_decision="HOLD_DATA_AUDIT; 不启动30结构或pilot"),
    audit_row(**common,
              system_id="MDMX_PMI_L10A_2021", priority="2 条件性同源受体复核但未放行",
              receptor_construct="human synthetic MDMX(24-108); 85 aa; 实际WT/AAA身份待定",
              construct_tags_or_mutations="5UML及3EQY相对O15151有Q68A/Q69A/E70A; 不静默当作assay构建体",
              affinity_evidence="SPR Kd nM: L10A387±22; T1L91±5; A4L53±4; N8A52±4; S11M106±3",
              enhancement_and_loss_evidence="增强 L10A→A4L/N8A; 减弱 T1L/S11M→L10A 为两个不同增强边反向; 每边仅计一次",
              PDB="5UML 3.0A; M3 peptide LTFLEYWAQLMQ; 3EQY相关PMI结构",
              structure_vs_assay_match="受体AAA与canonicalWT冲突未解释; 候选肽不是M3; 晶体肽链Q12坐标均缺失",
              training_overlap_audit="已查本地清单未见5UML/3EQY精确ID; 共享5UMM/PMI同源家族; 不等于训练无重叠",
              missing_information="实际SPR受体WT/AAA全序列; 肽及受体端化学; SPR误差定义; 新骨架独立证据",
              current_decision="HOLD_DATA_AUDIT; 不启动; 非第二无关家族"),
    audit_row(system_id="MDM2_MDMX_pDI_pDIQ_2010", priority="备用来源筛除",
              source_doi="10.1074/jbc.M109.073056", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2804373/",
              source_location="Methods ELISA/FP; Table1; local primary HTML",
              receptor_construct="ELISA GST-MDM2(1-150)或GST-MDMX(1-200); 结构构建体另有范围",
              construct_tags_or_mutations="GST assay fusion; pDI FP 使用荧光标记肽; 不与结构短受体等同",
              peptide_length="pDI/pDIQ 12 aa; pDIQ含四个替换; 其他突变逐条未进入输入",
              peptide_N_terminal="所选无标记突变面板未完成端化学核查; FP标记需另列",
              peptide_C_terminal="未完成逐条原始核查",
              assay_and_conditions="ELISA inhibition与FP不是同一观测; 不混用",
              quantity_reported="主要突变面板 IC50; 部分标记pDI的Kd不替代突变Kd",
              affinity_evidence="pDIQ ELISA IC50 MDM2 8nM/MDMX110nM; 不是Kd",
              error_and_replicates="不作为主Kd面板; 未为全部突变建立统一Kd误差口径",
              enhancement_and_loss_evidence="有增强抑制活性与减弱变化; 四突变pDIQ不是单点对照; 标签为IC50",
              PDB="论文有pDI/pDIQ晶体; 本轮未逐个下载/对齐; 不生成其结构输入",
              structure_vs_assay_match="GST大构建体与晶体短构建体不同; 不能直接拼接",
              used_in_current_project_development="本轮未评分; 完整逐变体曝光未审核; 不宣称盲测",
              training_overlap_audit="未作全面精确PDB清单审计; 未建立训练无重叠",
              selection_bias_or_independence="同旧MDM2家族且终点不符; 不补齐新Kd主面板",
              missing_information="同受体/同化学/同assay逐突变Kd及误差; 完整构建体/训练曝光审计",
              current_decision="EXCLUDE_PRIMARY_IC50_not_Kd"),
    audit_row(system_id="MCL1_BIM_2010", priority="新家族备用来源筛除",
              source_doi="10.1016/j.jmb.2010.03.058;10.1002/pro.329",
              source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2896288/;https://pmc.ncbi.nlm.nih.gov/articles/PMC2866276/",
              source_location="Dutta FP methods/competition and Fig3; Fire binding tables/methods; primary HTML",
              receptor_construct="Dutta human Mcl1(172-327), c-myc tag; Fire另有结合与晶体构建体; 禁止跨法合并",
              construct_tags_or_mutations="2PQK/3KZ0 receptor entity158aa含GS；受体原始结合试剂需逐法核对",
              peptide_length="Dutta synthetic Bim variants23aa; 2PQK entity29aa含GGSG及尾部; 3KZ0 MB7为23aa不同序列",
              peptide_N_terminal="Dutta大多数Ac-N; XG10/XD5为free-N例外; Fire各结合试剂未完成统一核查",
              peptide_C_terminal="Dutta大多数C-amide; XG10/XD5为free-acid例外; 不套用到Fire",
              assay_and_conditions="Dutta FP competition25C, FITC-Bim tracer; Fire直接FITC与competition构建体/方法另列",
              quantity_reported="Dutta Ki含删失; Fire有<2nM直接Kd及其他IC50; 不当作同assay点Kd",
              affinity_evidence="Dutta F4aV改善Ki<0.1nM相对Bim约2nM; 属于Ki界限而非精确Kd",
              error_and_replicates="删失界不能填为点值; 跨论文构建体/误差/标签不统一",
              enhancement_and_loss_evidence="SPOT与FP Ki有增强/减弱; 强亲和力Kd被检测下限删失; 本轮无合格同法Kd双向小面板",
              PDB="2PQK WT 2.0A;3KZ0 MB7 2.349A; Fire另有3KJ0/3KJ1/3KJ2",
              structure_vs_assay_match="已下载2PQK/3KZ0; resolved与entity不同; 结构不是F4aV候选或完全等长assay肽",
              used_in_current_project_development="本轮未评分; 新家族精确历史曝光未逐变体核查",
              training_overlap_audit="本地已查精确2PQK/3KZ0/3KJ0/1/2无命中; 同源训练重叠未排除",
              selection_bias_or_independence="可作将来新家族来源但终点/删失不满足冻结SPR-Kd主方案",
              missing_information="独立同构建体同assay可分辨的双向突变Kd与误差; 逐条标签/端化学",
              current_decision="EXCLUDE_PRIMARY_Ki_censoring_construct_mismatch"),
    audit_row(system_id="TIAM1_PDZ_SDC1", priority="新家族备用来源筛除",
              source_doi="10.3389/fmolb.2017.00065;10.1016/j.bpj.2018.01.008",
              source_url="https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2017.00065/full;https://pmc.ncbi.nlm.nih.gov/articles/PMC5883551/",
              source_location="Panel2017 Methods/tables; Panel2018 Table1/FEP provenance",
              receptor_construct="Tiam1 PDZ; WT与L920V等受体变体必须区分; 未释放exact输入",
              construct_tags_or_mutations="受体L920V增益不算肽单突变增强; 标签/边界逐原始实验尚未统一",
              peptide_length="Sdc1家族短末端肽; 6-8aa视具体实验/结构; 不套一个长度到全表",
              peptide_N_terminal="部分实验dansyl N标记; 不能与无标记链视为同化学",
              peptide_C_terminal="PDZ末端识别需free carboxyl; 逐试剂化学仍需原始实验确认",
              assay_and_conditions="FP/原始实验与PB-LIE/FEP预测区分; 不是表内全部均实测",
              quantity_reported="实测Kd/转化free energy与预测free energy混列; 不把A0V/A0mA预测标签当Kd",
              affinity_evidence="本轮未选择完整可比的实测Kd增强/减弱肽面板; 不编造数值",
              error_and_replicates="需回溯每条原始实验; simulation误差不等于Kd实验误差",
              enhancement_and_loss_evidence="WT-Sdc1大部分所查肽突变减弱; 受体增益和预测肽不能补实验增强",
              PDB="3KZE相关Sdc1;4GVD含标记化学; 未在本轮下载CIF",
              structure_vs_assay_match="未完成逐条terminal/length/tag/construct一致性审计",
              used_in_current_project_development="项目未作此来源新评分; 精确肽曝光未全面审核",
              training_overlap_audit="3kze_D本地Training_set/list line773; train.csv firstline30829 12rows",
              selection_bias_or_independence="非GraphPep训练无重叠; 不可把计算文献验证表全称实验独立数据",
              missing_information="同WT受体/同assay真实增强和减弱Kd肽; 端化学支持; 原始实验溯源",
              current_decision="EXCLUDE_PRIMARY_missing_experimental_gain_and_chemical_match"),
    audit_row(system_id="CAL_PDZ_iCAL36_2014", priority="新家族备用来源筛除",
              source_doi="10.1371/journal.pone.0103650", source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC4138078/",
              source_location="Table1/FP methods/results; primary local HTML title checked",
              receptor_construct="CAL PDZ CALP; exact boundary/tag逐assay未完整回溯",
              construct_tags_or_mutations="部分ITC用polyhistidine-tagged CALP; 不跨构建体合并",
              peptide_length="iCAL36 ANSRWPTSII 10aa及同长化学修饰类似物",
              peptide_N_terminal="native N与Ac-N明确区分; F* tracer含aminohexanoic-linker fluorescein",
              peptide_C_terminal="PDZ识别末端; 本轮未将全表逐条化学确定为输入",
              assay_and_conditions="FP竞争/ITC/SPR分别报告; 不能合并",
              quantity_reported="Table1表头Kd-from-FP; FP Methods实际说明Ki拟合; 不静默转换",
              affinity_evidence="FP iCAL36 22.6±8.0uM; Ac-N22.6±1.7; 多种Ac-K/halogen-K类似物; 不是标准AA单点小面板",
              error_and_replicates="FP三次separate stock experiments; 各方法误差需回溯; 本轮不设主Kd标签",
              enhancement_and_loss_evidence="有化学修饰增强/减弱; 部分增益小且非标准残基; 不满足当前三倍标准AA输入",
              PDB="论文列4NMO/4NMP/4NMQ/4NMR/4NMS/4NMT/4NMV; 本轮未下载对齐",
              structure_vs_assay_match="修改基团必须显式表达; 不转成普通K后继续评分",
              used_in_current_project_development="本轮未评分; 未完成该化学面板精确历史曝光核查",
              training_overlap_audit="未作全套本地训练PDB审计; 不宣称无重叠",
              selection_bias_or_independence="新家族但化学与测量定义不满足已冻结主方案",
              missing_information="精确construct和各试剂端化学/非标准残基支持; 统一同法实验Kd双向面板",
              current_decision="EXCLUDE_PRIMARY_modified_residue_and_quantity_definition"),
]
write_csv("candidate_source_audit.csv", rows)

source_rows: list[dict[str, object]] = []
for path in sorted((TASK / "sources").iterdir()):
    if not path.is_file():
        continue
    stem = path.stem
    title = doi = ""
    if path.suffix == ".html" and stem.startswith("PMC"):
        text = path.read_text(encoding="utf-8")
        for field, result_name in [("citation_title", "title"), ("citation_doi", "doi")]:
            match = re.search(r'<meta[^>]+name="' + field + r'"[^>]+content="([^"]+)"', text)
            if match:
                if result_name == "title": title = html.unescape(match.group(1))
                else: doi = match.group(1)
    if path.suffix == ".cif":
        url = f"https://files.rcsb.org/download/{stem}.cif"
        kind, status = "original_mmcif", "parsed_construct_only_no_scores"
    elif path.suffix == ".fasta":
        url = f"https://rest.uniprot.org/uniprotkb/{stem}.fasta"
        kind, status = "official_canonical_sequence", "construct_range_comparison"
    elif stem == "Frontiers_2017_00065":
        url = "https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2017.00065/full"
        doi = "10.3389/fmolb.2017.00065"
        kind, status = "original_paper_html", "reserve_source_not_primary_panel"
    elif path.suffix == ".txt":
        url = f"https://pmc.ncbi.nlm.nih.gov/articles/{stem}/"
        kind, status = "derived_plaintext_from_local_HTML", "navigation_only_verify_against_HTML"
    elif path.suffix == ".xml":
        accession = stem.split("_")[0]
        url = f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{accession}/unicode"
        kind, status = "BioC_retrieval_response", "not_used_for_primary_numeric_extraction"
        if accession == "PMC2866276":
            status = "not_validated_as_fulltext_retrieval_artifact_only"
    else:
        url = f"https://pmc.ncbi.nlm.nih.gov/articles/{stem}/"
        kind, status = "original_paper_html", "source_audit"
        if stem == "PMC8463443": status = "PRIMARY_TABLE2_SPR_NUMERIC_SOURCE"
    source_rows.append({"file": "sources/" + path.name, "source_url": url,
                        "source_kind": kind, "retrieval_date_local": "2026-10-07",
                        "paper_title": title, "doi": doi, "bytes": path.stat().st_size,
                        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "status": status})
write_csv("source_manifest.csv", source_rows)
print(f"Wrote {len(rows)} audited sources/systems, including two conditional candidates; {len(source_rows)} source artifacts.")
