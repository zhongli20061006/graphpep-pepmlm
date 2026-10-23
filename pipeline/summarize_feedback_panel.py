"""Create one row per peptide and descriptive robustness diagnostics."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/feedback_panel_mdm2_20261005"


def main():
    models = pd.read_csv(OUT / "analysis/model_results.csv")
    candidates = pd.read_csv(OUT / "candidates.csv", keep_default_na=False)
    seq = pd.read_csv(OUT / "analysis/sequence_summary.csv")
    metadata = json.loads((OUT / "parent_metadata.json").read_text())
    parents = {parent["peptide_id"]: parent for parent in metadata["parents"]}
    if len(models) != 117 or len(seq) != 36 or len(candidates) != 39:
        raise ValueError("Incomplete data for records and robustness diagnostics")
    records = []
    for candidate in candidates.to_dict("records"):
        group = models[models.peptide_id == candidate["peptide_id"]].sort_values("model_number")
        summary = seq[seq.peptide_id == candidate["peptide_id"]]
        is_parent = candidate["arm"] == "parent"
        parent = parents[candidate["parent_id"]]
        notes = []
        if parent["native_peptide_unresolved_positions"]:
            notes.append("parent crystal lacks last peptide residue; folding uses complete sequence")
        if len(group[group.structure_pass]) < 2:
            notes.append("insufficient quality-passing structures")
        records.append({"peptide_id": candidate["peptide_id"], "parent_id": candidate["parent_id"],
                        "sequence": candidate["sequence"], "length": candidate["length"], "arm": candidate["arm"],
                        "mask_positions": candidate["mask_positions"], "mutations": candidate["mutations"],
                        "PepMLM_ppl": candidate["ppl_score"], "parent_source_pdb": parent["pdb"],
                        "structure_prediction_method": "ColabFold 1.6.3 AFM-v3 + restrained Amber",
                        "afm_seed": metadata["seed"], "parameter_models": "1;2;3",
                        "ipTM_median": float(group.iptm.median()),
                        "peptide_pLDDT_mean": float(group.peptide_plddt.mean()),
                        "interface_PAE_median_A": float(group.interface_pae.median()),
                        "structure_pass_models": int(group.structure_pass.sum()),
                        "GraphPep_score_median_all_models": float(group.binding_score.median()),
                        "GraphPep_score_std_all_models": float(group.binding_score.std()),
                        "qualified_paired_gain_vs_parent": 0. if is_parent else summary.median_gain_vs_predicted_parent.iloc[0],
                        "predefined_improvement": False if is_parent else bool(summary.predefined_improvement.iloc[0]),
                        "residue_level_feedback_paths": ";".join(group.residue_feedback_path.dropna().astype(str)),
                        "complex_structure_paths": ";".join(group.complex_path.astype(str)),
                        "notes": "; ".join(notes)})
    pd.DataFrame(records).to_csv(OUT / "experiment_record.csv", index=False)
    stability = []
    for parent_id, group in models[models.arm != "parent"].groupby("parent_id"):
        valid = group.groupby("peptide_id").structure_pass.sum()
        complete = group[group.peptide_id.isin(valid[valid == 3].index)]
        matrix = complete.pivot(index="peptide_id", columns="model_number", values="binding_score")
        correlations = matrix.rank().corr()
        pairs = [correlations.loc[a, b] for a, b in [(1, 2), (1, 3), (2, 3)]]
        stability.append({"parent_id": parent_id, "children_all_3_models_passing": len(matrix),
                          "mean_pairwise_model_rank_correlation": float(np.nanmean(pairs)),
                          "min_pairwise_model_rank_correlation": float(np.nanmin(pairs)),
                          "median_child_score_std": float(matrix.std(axis=1).median())})
    pd.DataFrame(stability).to_csv(OUT / "analysis/score_stability.csv", index=False)
    seq["anchor_mutated"] = seq.anchor_positions_mutated.fillna("").astype(str).str.len() > 0
    diagnostic = seq.groupby(["parent_id", "arm", "anchor_mutated"]).agg(
        children=("peptide_id", "size"), improved=("predefined_improvement", "sum"),
        comparable_children=("median_gain_vs_predicted_parent", "count"),
        mean_qualified_gain=("median_gain_vs_predicted_parent", "mean")).reset_index()
    diagnostic.to_csv(OUT / "analysis/anchor_diagnostic_summary.csv", index=False)
    parent = parents["P53LIKE_3G03"]
    focus = seq[(seq.parent_id == "P53LIKE_3G03") & seq.arm.isin(["predicted_guided", "predicted_random"])]
    terminal_diagnostic = []
    for arm in ["predicted_guided", "predicted_random"]:
        group = focus[focus.arm == arm]
        s12_proline = group[group.mutations.str.contains("S12P", regex=False)]
        terminal_diagnostic.append({"parent_id": "P53LIKE_3G03", "arm": arm,
                                    "native_terminal_position_12_resolved": 12 not in parent["native_peptide_unresolved_positions"],
                                    "children_with_S12P": len(s12_proline),
                                    "S12P_children_meeting_improvement": int(s12_proline.predefined_improvement.sum())})
    pd.DataFrame(terminal_diagnostic).to_csv(OUT / "analysis/unresolved_terminal_mutation_check.csv", index=False)
    terminal_note = (f"P53LIKE 3G03中，预测反馈组有{terminal_diagnostic[0]['children_with_S12P']}条 S12P 子代，"
                     f"其中{terminal_diagnostic[0]['S12P_children_meeting_improvement']}条达标；对应随机组有"
                     f"{terminal_diagnostic[1]['children_with_S12P']}条 S12P 子代，其中"
                     f"{terminal_diagnostic[1]['S12P_children_meeting_improvement']}条达标。该父代实验结构未解析S12，"
                     "这个小样本观察不能证明预测反馈选位优于随机，也不能说明S12P的实测亲和力。")
    report_path = OUT / "report_zh.md"
    report = report_path.read_text(encoding="utf-8")
    report = report.split("\n## 逐肽记录与辅助诊断\n")[0]
    report += ("\n## 逐肽记录与辅助诊断\n\n"
               "`experiment_record.csv`有39行，每肽记录序列、长度、预测方法、置信度、GraphPep聚合分数、逐残基反馈路径和备注。`analysis/feedback_panel_gain.png`与SVG为改善分数图；每点是一条子代，数字为每组三条中达到预设标准的条数。\n\n"
               "`analysis/score_stability.csv`检查全部3个模型通过结构过滤的子代之间的评分排序一致性；`analysis/anchor_diagnostic_summary.csv`按是否改变F3/W7/L10分层。这两项为事后诊断，不能替代预设比较或消除小样本混杂。三条父代均属于同一MDM2体系，训练集排重尚未完成，本轮不构成未知靶点泛化验证。\n\n"
               f"事后检查：{terminal_note} 细节见`analysis/unresolved_terminal_mutation_check.csv`。\n\n"
               "GraphPep原论文以蛋白–肽docking/AlphaFold构象和接触准确性作为评价对象；本轮观察到的是评分代理上的变化。[GraphPep论文](https://www.nature.com/articles/s42256-025-01136-1)。\n")
    report_path.write_text(report, encoding="utf-8")
    print(pd.DataFrame(stability).to_string(index=False))
    print("Wrote 39 peptide records and descriptive robustness diagnostics")


if __name__ == "__main__":
    main()
