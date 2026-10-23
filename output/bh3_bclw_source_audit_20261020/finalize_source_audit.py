"""Finalize the source-only Foight 2014 BH3 audit; no prediction or scoring."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
CRITERIA = PROJECT / "output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(name: str, payload: object) -> None:
    (ROOT / name).write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    pairs = json.loads((ROOT / "paired_direction_candidates.json").read_text(encoding="utf-8"))
    overlap = json.loads((ROOT / "graphpep_overlap_audit.json").read_text(encoding="utf-8"))
    criteria_hash = sha256(CRITERIA)

    rows = []
    for pair in pairs["pairs"]:
        p = pair["parent_measurement"]
        v = pair["variant_measurement"]
        mutation = pair["mutation"]
        rows.append({
            "source": "Foight et al. 2014, DOI 10.1021/cb500340w, Table 1-2",
            "receptor": pair["receptor"],
            "parent": pair["parent"],
            "variant": pair["variant"],
            "mutation_local_1based": mutation["local_position_1based"],
            "mutation": f"{mutation['from']}{mutation['local_position_1based']}{mutation['to']}",
            "parent_Kd_nM_reported_95fitCI": p["raw"],
            "variant_Kd_nM_reported_95fitCI": v["raw"],
            "ratio_variant_over_parent": f"{pair['ratio_variant_over_parent']:.6g}",
            "point_direction": pair["point_effect_direction"],
            "approx_adverse_covariance_fit_interval_low": f"{pair['fit_conditional_approximate_log_ratio_interval_adverse_covariance'][0]:.6g}",
            "approx_adverse_covariance_fit_interval_high": f"{pair['fit_conditional_approximate_log_ratio_interval_adverse_covariance'][1]:.6g}",
            "marginal_endpoint_envelope_low_NOT_joint95": f"{pair['marginal_95CI_endpoint_envelope_NOT_joint95CI'][0]:.6g}",
            "marginal_endpoint_envelope_high_NOT_joint95": f"{pair['marginal_95CI_endpoint_envelope_NOT_joint95CI'][1]:.6g}",
            "fit_sensitivity_direction_supported": pair["direction_supported_under_declared_fit_sensitivity"],
            "formal_v2_label_released": False,
        })
    with (ROOT / "candidate_source_audit.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fit_v2 = {
        "version": "bh3_fit_summary_sensitivity_v2",
        "frozen_before_new_prediction_or_score": True,
        "source_semantics": "Direct Kd values in Foight Table 2; at least three replicate curves over three days, fit to averaged anisotropy data; plotted pointwise bars are SD; table uncertainty is 95% Matlab confint on nonlinear fit.",
        "ratio_method": "First-order delta-method interval on log(Kd_variant/Kd_parent), using reported 95% fit halfwidths as 1.96 local SEs. Report both independence and worst adverse covariance (sum of relative halfwidths) sensitivities. This is fit-conditional and approximate; it is not a source-reported ratio CI or a biological/replicate-level CI.",
        "printed_equation_audit": {
            "finding": "The printed direct-binding equation in Supplementary Methods appears to place 4*c*Kd under the square root, while the declared variables are peptide concentration c and receptor concentration x. The standard 1:1 mass-balance term is 4*c*x; as printed, the expression does not give zero bound fraction at x=0 for Kd>c.",
            "interpretation": "Likely a publication/typesetting error, but without raw anisotropy points or the Matlab fitting code we cannot establish whether the implemented fits used the printed or intended expression. Do not silently refit or treat this as proof that the published Kd values are wrong.",
            "visual_source": "figures/cb500340w_si_001_page4.png; inspected against the published supplement image",
        },
        "limits": [
            "No raw replicate-level anisotropy observations, per-run Kd fits, fit covariance, or Matlab code were found in the 16-page supplement.",
            "No division by sqrt(n); the three-day curve replicates are not three reported Kd estimates.",
            "The fitted peptide concentration was estimated from plateaued curves and reused; lower baseline was fixed.",
            "The ratio intervals omit uncertainty in active receptor concentration, oligomeric state, and fluorescence label effects.",
            "No formal benchmark labels are released by this sensitivity artifact; it preserves source candidates for pre-score review.",
        ],
    }
    write_json("fit_uncertainty_sensitivity_v2.json", fit_v2)

    raw_sources = [
        ("raw/PMC4168798.xml", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4168798/fullTextXML", "Europe PMC full-text XML; Foight et al. 2014"),
        ("raw/Foight2014_main.pdf", "https://dspace.mit.edu/server/api/core/bitstreams/9896bb9a-a783-4602-adc8-9ace8dfa944e/content", "MIT DSpace author/publication PDF"),
        ("raw/cb500340w_si_001.pdf", "https://ndownloader.figshare.com/files/3611775", "ACS Figshare supplement, article 2040459"),
        ("raw/PMC2896288_efetch.xml", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=2896288", "NCBI E-utilities XML; Dutta et al. 2010"),
        ("raw/Q92843_BCLW.fasta", "https://rest.uniprot.org/uniprotkb/Q92843.fasta", "UniProt canonical BCL-w sequence"),
        ("raw/Q07817_BCLXL.fasta", "https://rest.uniprot.org/uniprotkb/Q07817.fasta", "UniProt canonical BCL-xL sequence"),
    ]
    for pdb in ("1zy3", "3fdl", "4cim"):
        for kind in ("entry", "entity1", "entity2"):
            if kind == "entry":
                url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb}"
            else:
                url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{kind[-1]}"
            raw_sources.append((f"raw/rcsb_{pdb}_{kind}.json", url, f"RCSB PDB API {pdb.upper()} {kind}"))
    manifest_sources = []
    for rel, url, description in raw_sources:
        path = ROOT / rel
        manifest_sources.append({"path": rel, "url": url, "description": description, "bytes": path.stat().st_size, "sha256": sha256(path), "md5": md5(path) if rel.endswith("cb500340w_si_001.pdf") else None, "retrieval_status": "present_nonempty" if path.stat().st_size else "empty"})
    failed = []
    for rel, url, reason in [
        ("raw/figshare_si_response.html", "https://acs.figshare.com/articles/journal_contribution/2040459", "empty response artifact; authoritative Figshare API download was obtained separately"),
        ("raw/PMC4168798_supplementary.zip", "https://pmc.ncbi.nlm.nih.gov/articles/PMC4168798/", "empty/failed retrieval artifact; use verified Figshare PDF instead"),
    ]:
        path = ROOT / rel
        if path.exists():
            failed.append({"path": rel, "url": url, "bytes": path.stat().st_size, "status": "failed_or_empty_not_used_as_evidence", "note": reason})
    write_json("source_manifest.json", {
        "audit_date": "2026-10-20",
        "primary_sources": manifest_sources,
        "failed_or_empty_retrievals": failed,
        "figshare_supplement_expected_md5": "3d4b538a04191e436774e2d186977dfa",
        "figshare_supplement_md5_verified": md5(ROOT / "raw/cb500340w_si_001.pdf") == "3d4b538a04191e436774e2d186977dfa",
    })

    bclw = [p for p in pairs["pairs"] if p["receptor"] == "Bcl-w"]
    bclxl = [p for p in pairs["pairs"] if p["receptor"] == "Bcl-xL"]
    selected_gain = next(p for p in bclw if p["variant"] == "Bim_A2eT_I3dL")
    selected_loss = next(p for p in bclw if p["variant"] == "Bim_A2eT_I2dM")
    mapping = {
        "status": "DRAFT_PRE_SCORE_NOT_RUN",
        "scope": "Six-complex engineering pilot input representation only; not a full benchmark and not an assertion that model input exactly equals assay molecule.",
        "preferred_system": "Bcl-w fluorescence-anisotropy Kd, Foight et al. 2014 Table 2",
        "source_parent": "Bim_A2eT",
        "parent_sequence": "RPEIWITQELRRIGDEFNAYYAR",
        "gain_candidate": {"name": selected_gain["variant"], "sequence": "RPEIWITQELRRLGDEFNAYYAR", "mutation": "I3dL", "Kd_ratio_variant_over_parent": selected_gain["ratio_variant_over_parent"]},
        "loss_candidate": {"name": selected_loss["variant"], "sequence": "RPEIWMTQELRRIGDEFNAYYAR", "mutation": "I2dM", "Kd_ratio_variant_over_parent": selected_loss["ratio_variant_over_parent"]},
        "assay_target": {
            "peptide": "23-aa peptide; N-terminal 5/6-fluorescein amidite (isomer mixture as reported); C-terminal amide. Do not substitute the N-acetyl/C-amide competition peptides because those yield Ki, not the direct Kd labels.",
            "receptor": "human Bcl-w residues 1-164 with N-terminal c-myc tag after MBP/TEV purification; full c-myc/linker/TEV junction sequence is not recoverable from the cited methods and is left unknown.",
            "assay": "direct fluorescence anisotropy, 25 C, 25 mM Tris pH 7.8, 50 mM NaCl, 1 mM EDTA, 0.001% Triton X-100, 5% DMSO; 2 h read after equilibration check.",
        },
        "proposed_structure_model_input": {
            "receptor_chain": "UniProt Q92843 residues 1-164 only; omit c-myc/linker because its exact sequence is unavailable and it is a common terminal tag, while preserving this as an explicit assay-to-model mismatch.",
            "peptide_chain": "standard 23-aa sequence only; record N-fluorescein and C-amide in metadata, but do not claim those chemical groups were structurally modeled.",
            "interpretation": "Any structural scores are proxies evaluated against the labeled-peptide/c-myc-receptor assay Kd direction; they are not predictions of unlabeled-peptide Kd or complete binding free energy.",
        },
        "seed_policy": {"seeds": [20261008, 20261009], "complexes_max": 6, "purpose": "engineering feasibility and cross-seed stability only"},
        "entry_caveat": "Do not activate until this common-tag/core-only mapping is explicitly frozen in the run configuration and the source equation discrepancy is retained in the analysis record.",
    }
    write_json("input_representation_draft.json", mapping)
    frozen_mapping = {
        **mapping,
        "status": "FROZEN_PRE_SCORE_ENGINEERING_PILOT_MAPPING",
        "entry_caveat": "Mapping frozen for the limited engineering pilot before any new score. This is an explicit core-only proxy of the labeled assay system; it does not make the modeled molecule chemically identical to the assay molecule. Keep this mismatch attached to every result.",
        "score_or_prediction_seen_at_freeze": False,
    }
    write_json("input_representation_frozen_pilot_v1.json", frozen_mapping)

    counts = {}
    for receptor, subset in (("Bcl-w", bclw), ("Bcl-xL", bclxl)):
        counts[receptor] = {"source_defined_single_substitution_pairs": len(subset), "gain_point_candidates": sum(x["point_effect_direction"] == "gain" for x in subset), "loss_point_candidates": sum(x["point_effect_direction"] == "loss" for x in subset), "fit_sensitivity_supported": sum(bool(x["direction_supported_under_declared_fit_sensitivity"]) for x in subset)}
    decision = {
        "version": "foight_bh3_source_decision_v1",
        "frozen_v2_criteria_sha256": criteria_hash,
        "new_prediction_or_score_seen": False,
        "source_screen": "QUALIFIES_AS_BIDIRECTIONAL_SOURCE_CANDIDATE; public and development-exposed, not blind",
        "preferred_small_pilot_system": "Bcl-w",
        "pilot_recommendation": "CONDITIONAL_GO_TO_RUN_CONFIGURATION_FREEZE_FOR_SIX_COMPLEX_ENGINEERING_PILOT",
        "pilot_input_mapping_status": "FROZEN_PRE_SCORE_AS_AN_EXPLICIT_CORE_ONLY_PROXY; model input is not chemically identical to the labeled assay molecule",
        "pilot_scope": "Bim_A2eT parent, I3dL gain candidate, I2dM loss candidate; seeds 20261008 and 20261009; six complexes maximum; no method-efficacy claim.",
        "full_panel_decision": "NO_GO_FROM_THIS_PAPER_ALONE",
        "systems": counts,
        "why_bclw_first": [
            "Same Bim_A2eT parent, same receptor and direct Kd assay; source-defined gain I3dL (12 vs 39 nM; ratio 0.308) and loss I2dM (260 vs 39 nM; ratio 6.67).",
            "The candidate mutation directions remain separated from 1 under the declared fit-conditional adverse-covariance sensitivity, while the fit limitations and printed-equation discrepancy remain explicit.",
            "GraphPep's packaged structures include Bcl-w-family exposure (4CIM, mutated receptor/other BH3 partner) but the audit found no exact Bim_A2eT mutant peptide; this is structural exposure, not evidence that Table 2 Kd labels trained GraphPep.",
        ],
        "why_not_full_benchmark": [
            "Each receptor contributes only seven source-defined pairs; Bcl-w has two gains and five losses, Bcl-xL one gain and six losses. Neither reaches >=10 pairs and >=4 per class, and pooling receptors is prohibited.",
            "GraphPep packaged training includes the Bcl-xL:Bim 4QVF complex; model-weight/pretraining sources are not exhaustively known. Do not call this a blind test.",
        ],
        "material_limits": [
            "Paper Table 2 is direct Kd, but reports a fit confidence interval from averaged anisotropy curves rather than independent Kd estimates; no raw curves or MATLAB code were found.",
            "The supplement's printed quadratic binding equation appears inconsistent with the declared variables; likely typographical, but not independently resolvable from published raw data.",
            "The complete c-myc/TEV junction is unspecified; the assay reports Bcl-w 1-164 plus N-terminal c-myc tag. The proposed model input omits that common tag and the fluorescein label and is explicitly a proxy.",
            "Dutta et al. 2010 reported substantial Bcl-w oligomer in their preparation; this is a protocol-level warning, not a measured oligomer fraction for the Foight 2014 batch.",
        ],
        "not_done": ["no new structure prediction", "no scorer execution", "no model training", "no paid service or external contact", "no wet-lab order"],
    }
    write_json("decision.json", decision)

    protocol = f'''# Foight 2014 BH3：六结构工程试跑应用稿

状态：**输入映射已在新分数出现前冻结；预测与评分尚未启动。** 这是冻结判据的来源应用稿，不修改`direction_criteria_v2`。唯一对应的机器可读文件是`input_representation_frozen_pilot_v1.json`。

## 目标与候选

只检查复合物预测是否可重复，以及两枚固定seed是否造成方向翻转；不据此宣称任何评分器能预测真实亲和力。主体系选Bcl-w直接Kd实验：亲本为工程序列`Bim_A2eT`（`RPEIWITQELRRIGDEFNAYYAR`）；增益候选`I3dL`（`RPEIWITQELRRLGDEFNAYYAR`，Kd比值{selected_gain['ratio_variant_over_parent']:.3f}）；减弱候选`I2dM`（`RPEIWMTQELRRIGDEFNAYYAR`，Kd比值{selected_loss['ratio_variant_over_parent']:.3f}）。所有比较都回到同一亲本，不跨受体汇总。

原文表2报告亲本39±9.4 nM、I3dL 12±2.6 nM、I2dM 260±52 nM；±为95%曲线拟合区间，不是SD/SEM。作者将至少三天、至少三次曲线先平均后拟合。局部不确定度文件仅作拟合条件近似敏感性分析，不等于原文提供了突变比值95%区间。

## 输入映射

- assay目标：23 aa肽带N端5/6-fluorescein amidite、C端酰胺；受体为人Bcl-w 1–164加N端c-myc标签，MBP/TEV纯化后测量。直接荧光各向异性Kd，不拿N-acetyl竞争实验的Ki替代。
- 建模输入：Bcl-w核心1–164和标准氨基酸肽序列。模型不编码共同的c-myc/TEV接头、荧光染料或酰胺端；这项差异必须保留在结果表，结论只针对该proxy与论文Kd方向的一致性。
- 补充材料把结合方程印成含`4*c*Kd`的根号项，按其变量定义通常应为`4*c*x`。将此列为文献复现风险，不重拟合、不修正原值。

## 运行边界

最大6个复合物：3条肽×AFM seed 20261008、20261009；模型1、3 recycles、关闭templates、共用受体MSA。预算上限2 GPU小时、2 GiB；无需PepMLM生成、训练或全池运行。首轮仅评估结构生成成功率、链/残基映射、接触是否存在、seed间pose及预冻结评分差方向稳定性。

FoldX、GraphPep、ESM2和接触数保持同等候选地位，必须先记录各自精确版本/权重哈希/命令/聚合和失败规则；低置信、无接触、不可评分均保留为失败。GraphPep越低越好；parent−variant为正只表示代理分数改善。该pilot不能证明reward有效，也不能替代≥10对、增益/减弱各≥4且跨至少两个亲本背景的完整面板。

## 停止条件

输入映射未冻结、亲本/突变链映射不一致、2小时或2 GiB预算将超限时停止。失败复合物不删去、不改seed、不重抽亲本。即使两seed稳定，也只报告工程可行性；等来源面板满足冻结门后再计算主指标。
'''
    (ROOT / "pilot_protocol_application_zh.md").write_text(protocol, encoding="utf-8")

    report = f'''# Foight 2014 BH3来源审计

日期：2026-10-20。**资料与输入审计完成；未运行新结构预测、模型评分或训练。** 这篇论文给出了独立于MDM2/PMI开发选择的可用突变方向候选，因此主线继续推进；它没有提供足够样本组成完整验证面板，也不能称盲测。

## 结论先行

建议将**Bcl-w / Bim_A2eT**作为六结构工程试跑的首选来源：同一亲本下，I3dL为增强候选、I2dM为减弱候选；报告点估计比分别为{selected_gain['ratio_variant_over_parent']:.3f}和{selected_loss['ratio_variant_over_parent']:.2f}。这是来源文献中的真实Kd方向证据，不是模型效果证据。新分数出现前已冻结[输入表示](input_representation_frozen_pilot_v1.json)：模型使用“Bcl-w核心+标准肽序列”，不编码实验中的荧光标记、酰胺端和c-myc标签；所有分数都必须解释为对带标记assay Kd的结构代理。

**值得进入六结构工程pilot的运行配置冻结；不值得把这篇单独扩成完整方向benchmark。** 运行前仍须完成v2要求的命令、版本、模型/权重哈希和输出路径冻结。此处没有把该建议等同于已运行。

## 候选亲和力数据

Foight等在同一荧光各向异性直接结合assay中测量Bim家族23-mer肽对多个Bcl-2受体的Kd。`Bim_A2eT`是工程亲本，不是天然WT。Table 1明示变体关系；Table 2提供同亲本、同受体的定量Kd。对Bcl-w，亲本39±9.4 nM；I3dL 12±2.6 nM（变体/亲本=0.308，约3.25倍增强）；I2dM 260±52 nM（比值6.67，减弱）。

每个±是Matlab非线性曲线拟合的95%置信区间。补充材料说直接结合曲线由至少三天、至少三次重复先取平均再拟合，图中误差棒是SD；这些不是三份独立Kd，也不能除以√3。我们预先按源CI半宽作了拟合条件的log比值近似，另用最不利协方差半径作敏感性检查；这不是原文给出的突变比值95%CI，也不覆盖活性受体浓度、聚集或荧光标签影响。数据表把近似与边际端点包络分开保存，正式v2标签暂未释放。

| 受体 | 同背景源定义单点对 | 点估计增强/减弱 | pilot方向候选 |
|---|---:|---:|---|
| Bcl-w | 7 | 2 / 5 | Bim_A2eT → I3dL（增强）；→ I2dM（减弱） |
| Bcl-xL | 7 | 1 / 6 | 同亲本I3dL增强（比值0.447）；F4aI减弱（比值6.47） |

两受体不能合并凑样本。每个只有7对；都未达到完整面板“至少10对、每类至少4对、至少两个亲本背景且两类信号不由单一背景贡献”的门。Bcl-xL可作次选，但它的增强更靠近不确定区间边界，且GraphPep训练集包含4QVF Bcl-xL:Bim结构。

## 构建体、化学与assay

直接Kd用Table 1所列23-mer，N端5/6-fluorescein amidite，C端酰胺。论文另有未标记、N-acetyl/C-amide肽竞争试验，报告的是Ki；不与Table 2直接Kd混用。受体按Dutta等方法表达：Bcl-w 1–164、Bcl-xL 1–209，MBP融合经TEV切割后带N端c-myc标签。完整c-myc/TEV接头序列没有在这两篇方法中列出，所以当前结构输入采用核心序列代理，并明确省略该共同标签。

论文以2小时各向异性数据为主，并用24小时曲线检查平衡。缓冲液含25 mM Tris pH 7.8、50 mM NaCl、1 mM EDTA、0.001% Triton X-100和5% DMSO，温度25°C。论文指出荧光基团会影响部分肽的亲和力；所以这些数值应解释为带标记肽在该assay中的Kd，不等于未标记肽Kd。

Dutta 2010在其制备中报告Bcl-w有明显寡聚体比例。该数值不应挪作Foight 2014这批蛋白的测量结果，但提示用单体核心结构解释Bcl-w Kd时有构象/活性分数不确定性。它降低了Bcl-w作为完整物理亲和力benchmark的纯净度，却不抹去同一批assay内观察到的配对方向。

## 结构与开发暴露

- 1ZY3是Bcl-w与BID 20-mer的NMR/建模结构；受体长度/尾端构建和肽都不等于本论文Bcl-w/Bim样品。它只能作家族结构参考。
- 3FDL是Bcl-xL与Bim 26-mer结构，受体为截短/工程构建；不是论文中1–209+c-myc与荧光23-mer的精确assay输入。
- GraphPep随附5018份原生结构文件全部完成哈希核对：无缺失、无哈希偏差；初筛命中7条受体链、7条天然Bim肽链。Bcl-xL:Bim的4QVF复合物在GraphPep训练拆分；Bcl-w的4CIM训练结构是突变Bcl-w与另一BH3肽。此为结构/序列开发暴露，不证明Table 2的Kd标签用于GraphPep训练。Bim天然肽也出现在训练和验证结构中； engineered Bim_A2eT及所选四个子代未在该快照原生肽序列中命中。权重的全部预训练来源仍未知。
- Foight论文公开发表并详细讨论这些突变；本项目已经查看标签后选择来源，因此绝不称盲测。可称“独立于MDM2/PMI选择的公开跨家族外部文献面板”，并单独披露结构训练重叠。

## 复现风险与处理

补充材料的直接结合方程在印刷图像中似乎写为根号项`(c+x+Kd)^2−4*c*Kd`，而文字把`c`定义成肽浓度、`x`定义成受体浓度；通常的1:1质量平衡项应为`4*c*x`。按印刷式在x=0时可能给出非零结合。最可能是排版错误，但原始曲线和Matlab代码未公开，无法确认论文实际计算所用公式。我们保留这一复现风险，不改写、不重拟合、不判定表中Kd错误。

## 下一步建议

当前判断为**值得进入运行配置冻结**：Bcl-w亲本/I3dL/I2dM×两个既定seed，最多6个复合物、2 GPU小时、2 GiB。完整AFM命令、环境、权重/输出哈希和各候选scorer配置仍须在运行前冻结；本次只审来源与冻结输入映射，没有启动预测。先验证链映射和跨seed结构/分数可用性，再决定是否投资寻找更多合格突变数据。即使pilot稳定，也只证明工程输入可重复；不能说明GraphPep/FoldX/ESM的突变方向准确，更不等于Kd、ΔG(binding)或新设计肽亲和力得到验证。

原始主文、补充PDF、Dutta全文XML、RCSB/UniProt记录、逐项候选数据、GraphPep重叠清单及哈希见本目录`raw/`、`candidate_source_audit.csv`、`graphpep_overlap_audit.json`和`source_manifest.json`。冻结v2原判据SHA-256为`{criteria_hash}`，未改动。机器可读决定见`decision.json`。
'''
    (ROOT / "report_zh.md").write_text(report, encoding="utf-8")

    verification = {
        "source_table2_finite_and_censored_cells": json.loads((ROOT / "raw_table2_measurements.json").read_text(encoding="utf-8"))["cell_count"],
        "source_sequences": len(json.loads((ROOT / "peptide_sequences_source.json").read_text(encoding="utf-8"))["records"]),
        "source_defined_single_mutation_edges": len(json.loads((ROOT / "peptide_sequences_source.json").read_text(encoding="utf-8"))["edges"]),
        "paired_comparisons": len(pairs["pairs"]),
        "separate_system_counts": counts,
        "graphpep_native_files_examined": overlap["native_files_examined"],
        "graphpep_missing_native_files": len(overlap["missing_paths"]),
        "graphpep_manifest_hash_mismatches": len(overlap["manifest_sha256_mismatches"]),
        "graphpep_recB_files_included": overlap["native_path_name_counts"]["recb.pdb"],
        "criteria_v2_sha256_current": criteria_hash,
        "criteria_v2_unchanged_from_frozen": criteria_hash == "23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b",
        "figshare_supplement_md5_matches_api": md5(ROOT / "raw/cb500340w_si_001.pdf") == "3d4b538a04191e436774e2d186977dfa",
        "new_predictions_or_scores_seen": False,
    }
    write_json("verification.json", verification)

    excluded = {"artifact_manifest_sha256.json"}
    artifacts = []
    for path in sorted(ROOT.rglob("*")):
        if path.is_file() and path.name not in excluded:
            artifacts.append({"path": path.relative_to(ROOT).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_json("artifact_manifest_sha256.json", {"scope": "all files in this audit directory except this self-referential manifest", "artifacts": artifacts})
    print(json.dumps({"counts": counts, "verified": verification, "artifacts_hashed": len(artifacts)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
