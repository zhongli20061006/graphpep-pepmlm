"""Re-evaluate the previous seed under the panel's uniform preprocessing."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import feedback_panel as panel

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / "output/feedback_pilot_mdm2_20261005"
OUT = ROOT / "output/feedback_seed_retest_mdm2_20261005"


def score() -> None:
    metadata = json.loads((panel.OUT / "parent_metadata.json").read_text())
    parent = next(p for p in metadata["parents"] if p["peptide_id"] == "P53WT_1YCR")
    old_protocol = json.loads((PRIOR / "protocol.json").read_text())
    if panel.digest(PRIOR / "candidates.csv") != old_protocol["candidate_file_sha256"]:
        raise ValueError("Previous seed's frozen candidate sequences changed")
    if panel.digest(PRIOR / "msa_inputs/PARENT.a3m") != panel.digest(panel.OUT / "parent_msa/P53WT_1YCR.a3m"):
        raise ValueError("Parent receptor/peptide MSA differs between seeds")
    candidate = pd.read_csv(PRIOR / "candidates.csv", keep_default_na=False)
    specification = {**parent, "gate": panel.COMPARISON_PROTOCOL["model_gate"]}
    native_contacts = panel.collect_native_contact_positions(parent)
    model_rows = []
    for row in candidate.to_dict("records"):
        original_id = row["peptide_id"]
        row.update(parent_id="P53WT_1YCR", peptide_id=f"OLDSEED_{original_id}",
                   arm={"guided": "native_guided", "random": "native_random", "parent": "parent"}[row["arm"]])
        paths = sorted((OUT / "relaxed").glob(f"{original_id}_relaxed_rank_*_model_*_seed_20261005.pdb"))
        if len(paths) != 3:
            raise ValueError(f"Missing previous-seed relaxed structures: {original_id}")
        if row["arm"] != "parent":
            panel.common.validate_child(parent["sequence"], row["sequence"], [int(p) for p in str(row["mask_positions"]).split(";")])
        for path in paths:
            result = panel.metrics(specification, row, path, Path(parent["native_path"]), native_contacts)
            result["complex_path"] = result.pop("model_path")
            model_rows.append(result)
    frame = panel.run_graphpep(pd.DataFrame(model_rows), OUT / "graphpep")
    if len(frame) != 21:
        raise ValueError("Expected 21 previous-seed evaluator outputs")
    for row in frame.to_dict("records"):
        panel.score_bounds(row)
    frame.to_csv(OUT / "model_results.csv", index=False)
    baseline = frame[frame.arm == "parent"].set_index("model_number")
    summaries = []
    for ident, group in frame[frame.arm != "parent"].groupby("peptide_id", sort=False):
        group = group.set_index("model_number").sort_index()
        reference = baseline.loc[group.index]
        usable = group.structure_pass.astype(bool) & reference.structure_pass.astype(bool)
        child, reference = group[usable], reference[usable]
        gain = reference.binding_score - child.binding_score
        median_gain = float(reference.binding_score.median() - child.binding_score.median()) if len(child) >= 2 else np.nan
        improved = bool(len(child) >= 2 and median_gain > .05 and int((gain > 0).sum()) >= 2)
        summaries.append({"peptide_id": ident, "arm": group.arm.iloc[0], "sequence": group.sequence.iloc[0],
                          "mask_positions": group.mask_positions.iloc[0], "mutations": group.mutations.iloc[0],
                          "graphpep_median": float(group.binding_score.median()),
                          "quality_pass_models": int(group.structure_pass.sum()),
                          "common_passing_parent_child_models": len(child),
                          "median_gain_vs_predicted_parent": median_gain,
                          "predefined_improvement": improved, "matched_models_improved": int((gain > 0).sum()),
                          "anchor_positions_mutated": ";".join(str(p) for p in [3, 7, 10] if p in [int(s) for s in str(group.mask_positions.iloc[0]).split(";")])})
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "sequence_summary.csv", index=False)
    arms = summary.groupby("arm", sort=False).agg(n_children=("peptide_id", "size"),
             improved=("predefined_improvement", "sum"), mean_qualified_gain=("median_gain_vs_predicted_parent", "mean")).reset_index()
    arms.to_csv(OUT / "arm_summary.csv", index=False)
    panel.write_json(OUT / "verification.json", {
        "seed": 20261005, "structure_count": 21, "evaluator_rows": 21,
        "candidate_sequences_sha256": panel.digest(PRIOR / "candidates.csv"),
        "identical_parent_msa_between_seeds": True,
        "all_chain_sequences_and_graphpep_formulas_verified": True,
        "parent_passing_models": int(baseline.structure_pass.sum()),
        "preprocessing": panel.COMPARISON_PROTOCOL["relaxation"],
    })
    print(arms.to_string(index=False), flush=True)


def compare() -> None:
    old = pd.read_csv(OUT / "sequence_summary.csv")
    new = pd.read_csv(panel.OUT / "analysis/sequence_summary.csv")
    new = new[(new.parent_id == "P53WT_1YCR") & new.arm.isin(["native_guided", "native_random"])]
    combined = pd.concat([old.assign(seed=20261005), new.assign(seed=20261009)], ignore_index=True)
    combined.to_csv(OUT / "cross_seed_sequences.csv", index=False)
    arm = combined.groupby(["seed", "arm"]).agg(n_children=("peptide_id", "size"),
            improved=("predefined_improvement", "sum"), common_model_comparable=("median_gain_vs_predicted_parent", "count"),
            mean_qualified_gain=("median_gain_vs_predicted_parent", "mean")).reset_index()
    arm.to_csv(OUT / "cross_seed_arm_summary.csv", index=False)
    overlap = sorted(set(old[old.arm == "native_guided"].sequence) & set(new[new.arm == "native_guided"].sequence))
    lines = ["# p53父代跨种子复验（相同Amber预处理）", "",
             "两轮使用相同受体/父代MSA、AFM-v3参数模型1/2/3、3次recycle、相同结构过滤和GraphPep改善门槛。旧轮使用其冻结的原有6条子代，新轮使用新的采样种子。所有结构使用相同Amber松弛，旧轮的原始报告与评分保留。", "",
             "|AFM种子|选位策略|子代|通过改善标准|可用配对模型≥2的子代|平均合格配对分数增益|",
             "|---:|---|---:|---:|---:|---:|"]
    for row in arm.to_dict("records"):
        lines.append(f"|{row['seed']}|{row['arm']}|{row['n_children']}|{row['improved']}|{row['common_model_comparable']}|{row['mean_qualified_gain']:+.4f}|")
    lines += ["", f"两轮native-guided重复序列：{', '.join(overlap) if overlap else '无'}。若重复，只能视为换AFM种子的结构重验，不能算独立序列重复。", "",
              "每组仍仅3条子代；3个AFM参数模型不是独立生物学重复，GraphPep改善不等同于实验亲和力提高。此复验还同时改变PepMLM采样与AFM随机种子，因此不能把差异单独归因到某一个随机来源。", "",
              "各序列见`cross_seed_sequences.csv`；全部结构评价见`model_results.csv`；输入和公式核验见`verification.json`。"]
    (OUT / "report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(arm.to_string(index=False), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["score", "compare"])
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    score() if args.stage == "score" else compare()
