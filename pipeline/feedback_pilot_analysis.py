#!/usr/bin/env python3
"""Evaluate all pilot structures using frozen gates and the installed GraphPep."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer

import feedback_pilot as pilot
import step3_graphpep as gp

OUT = pilot.OUT
ANALYSIS = OUT / "analysis"
PARSER = PDBParser(QUIET=True)


def chain_data(structure, chain: str) -> tuple[str, list, np.ndarray, np.ndarray]:
    residues = [r for r in structure[0][chain] if r.resname in gp.THREE_TO_ONE]
    sequence = "".join(gp.THREE_TO_ONE[r.resname] for r in residues)
    xyz, positions = [], []
    for pos, residue in enumerate(residues, 1):
        for atom in residue:
            if atom.element not in {"H", "D"}:
                xyz.append(atom.coord)
                positions.append(pos)
    return sequence, residues, np.asarray(xyz, dtype=float), np.asarray(positions)


def contacts(rec_data: tuple, pep_data: tuple, cutoff: float) -> tuple[set, float, int]:
    xyz_r, pos_r = rec_data[2:]
    xyz_p, pos_p = pep_data[2:]
    distance = np.sqrt(np.sum((xyz_r[:, None, :] - xyz_p[None, :, :]) ** 2, axis=2))
    a, b = np.where(distance <= cutoff)
    pairs = {(int(pos_r[i]), int(pos_p[j])) for i, j in zip(a, b)}
    return pairs, float(distance.min()), int((distance < 1.5).sum())


def collect_metrics(config: dict, candidates: pd.DataFrame) -> pd.DataFrame:
    native = PARSER.get_structure("native", str(pilot.PRIOR / "POSCTRL_MDM2_P53_native.pdb"))
    native_rec, native_pep = chain_data(native, "A"), chain_data(native, "B")
    native_pairs, _, _ = contacts(native_rec, native_pep, 5.)
    nr, np_ = len(config["receptor_sequence"]), len(config["parent_sequence"])
    rows, sets = [], {}
    for candidate in candidates.to_dict("records"):
        ident = candidate["peptide_id"]
        files = sorted((OUT / "structures").glob(f"{ident}_unrelaxed_rank_*.pdb"))
        if len(files) != 3:
            raise ValueError(f"Expected 3 structures for {ident}, found {len(files)}")
        observed_models = set()
        for path in files:
            match = re.search(r"_model_(\d+)_seed_(\d+)\.pdb$", path.name)
            if not match:
                raise ValueError(f"Invalid AFM filename: {path.name}")
            model, seed = map(int, match.groups())
            if seed != config["afm_seed"] or model in observed_models:
                raise ValueError("Unexpected seed or repeated model")
            observed_models.add(model)
            structure = PARSER.get_structure(ident, str(path))
            rec, pep = chain_data(structure, "A"), chain_data(structure, "B")
            if rec[0] != config["receptor_sequence"] or pep[0] != candidate["sequence"]:
                raise ValueError(f"Chain sequence mismatch: {path}")
            if [c.id for c in structure[0]] != ["A", "B"]:
                raise ValueError("Unexpected chain order or additional chains")
            conf_path = path.with_name(path.name.replace("_unrelaxed_", "_scores_")).with_suffix(".json")
            conf = json.loads(conf_path.read_text())
            pae, plddt = np.asarray(conf["pae"]), np.asarray(conf["plddt"])
            if pae.shape != (nr + np_, nr + np_) or plddt.shape != (nr + np_,):
                raise ValueError("Confidence dimensions do not match chain sequences")
            if not np.isfinite(pae).all() or not np.isfinite(plddt).all():
                raise ValueError("Nonfinite AFM confidence values")
            pairs5, minimum, clashes = contacts(rec, pep, 5.)
            pairs6, _, _ = contacts(rec, pep, 6.)
            contacting_rec = sorted({r - 1 for r, _ in pairs6})
            pindex = list(range(nr, nr + np_))
            interface_pae = float(np.median(np.concatenate([
                pae[np.ix_(contacting_rec, pindex)].ravel(), pae[np.ix_(pindex, contacting_rec)].ravel()
            ]))) if contacting_rec else float("nan")
            sup = Superimposer()
            sup.set_atoms([r["CA"] for r in native_rec[1]], [r["CA"] for r in rec[1]])
            rotation, translation = sup.rotran
            aligned = np.asarray([r["CA"].coord for r in pep[1]]) @ rotation + translation
            native_ca = np.asarray([r["CA"].coord for r in native_pep[1]])
            peptide_rmsd = float(np.sqrt(np.mean(np.sum((aligned - native_ca) ** 2, axis=1))))
            gates = config["model_gate"]
            reasons = []
            if float(conf["iptm"]) < gates["iptm_min"]:
                reasons.append("low_iptm")
            pep_plddt = float(plddt[nr:].mean())
            if pep_plddt < gates["peptide_plddt_min"]:
                reasons.append("low_peptide_plddt")
            if not math.isfinite(interface_pae) or interface_pae > gates["interface_pae_max"]:
                reasons.append("high_or_missing_interface_pae")
            if minimum < gates["minimum_interchain_heavyatom_distance_min"]:
                reasons.append("heavyatom_clash")
            if len(pairs6) < gates["minimum_contacts"]:
                reasons.append("no_contact")
            key = (ident, model)
            sets[key] = pairs5
            rows.append({**candidate, "peptide_sequence": candidate["sequence"], "model_id": f"model_{model}",
                         "model_number": model, "afm_seed": seed, "complex_path": str(path),
                         "receptor_sequence": config["receptor_sequence"], "iptm": float(conf["iptm"]),
                         "ptm": float(conf["ptm"]), "peptide_plddt": pep_plddt,
                         "interface_pae": interface_pae, "minimum_heavyatom_distance": minimum,
                         "clashing_atompairs_under_1p5A": clashes, "contact_pairs_5A": len(pairs5),
                         "contact_pairs_6A": len(pairs6), "native_contact_recall_5A": len(pairs5 & native_pairs) / len(native_pairs),
                         "peptide_CA_RMSD_after_receptor_alignment_A": peptide_rmsd,
                         "structure_pass": not reasons, "quality_fail_reasons": ";".join(reasons),
                         "structure_sha256": pilot.sha256(path), "confidence_sha256": pilot.sha256(conf_path),
                         "confidence_path": str(conf_path)})
        if observed_models != set(config["afm_models"]):
            raise ValueError("Incorrect AFM model set")
    for row in rows:
        current = sets[(row["peptide_id"], row["model_number"])]
        overlaps = []
        for model in config["afm_models"]:
            if model == row["model_number"]:
                continue
            other = sets[(row["peptide_id"], model)]
            union = current | other
            overlaps.append(len(current & other) / len(union) if union else 0.)
        row["contact_jaccard_other_models"] = float(np.mean(overlaps))
    return pd.DataFrame(rows).sort_values(["peptide_id", "model_number"]).reset_index(drop=True)


def verify_graphpep_score(row: dict) -> dict:
    """Bound the formula using rounding intervals of the exported contact probabilities."""
    if row["graphpep_status"] == "no_contact":
        if int(row["contact_pairs_6A"]) != 0:
            raise ValueError("Unexpected missing contact score")
        return {"formula_verified": False, "formula_note": "no_contact_unscored"}
    values = pd.read_csv(row["contact_feedback_path"])["CONFIDENT_SCORES"].to_numpy()
    lower_count = sum(weight * int((values - .00005 > threshold).sum()) for weight, threshold in [(1, .5), (2, .7), (4, .9)]) / 7
    upper_count = sum(weight * int((values + .00005 > threshold).sum()) for weight, threshold in [(1, .5), (2, .7), (4, .9)]) / 7
    low_score, high_score = -math.log1p(upper_count), -math.log1p(lower_count)
    if not low_score - .000051 <= row["binding_score"] <= high_score + .000051:
        raise ValueError(f"GraphPep formula mismatch: {row['peptide_id']} {row['model_id']}")
    if len(values) != row["contact_pairs_6A"] or len(values) != row["n_contacts"]:
        raise ValueError("GraphPep contacts disagree with independent geometry")
    return {"formula_verified": True, "formula_score_lower_bound": low_score, "formula_score_upper_bound": high_score,
            "rounded_contacts_at_threshold": int(np.isin(values, [.5, .7, .9]).sum())}


def summarize(frame: pd.DataFrame, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    parent = frame[frame.arm == "parent"].set_index("model_number")
    if len(parent) != 3:
        raise ValueError("Missing predicted parent baseline")
    parent_median = float(parent.binding_score.median())
    parent_pass = int(parent.structure_pass.sum()) >= 2
    rows = []
    for ident, group in frame.groupby("peptide_id", sort=False):
        aligned_group = group.set_index("model_number").reindex(parent.index)
        gain = parent.binding_score - aligned_group.binding_score
        quality = group[group.structure_pass]
        median = float(group.binding_score.median())
        sequence_pass = len(quality) >= 2
        improvements = int((gain > 0).sum())
        improved = bool(parent_pass and sequence_pass and parent_median - median > config["minimum_score_gain"] and improvements >= 2)
        rows.append({"peptide_id": ident, "arm": group.arm.iloc[0], "sequence": group.sequence.iloc[0],
                     "mutations": group.mutations.iloc[0], "ppl_score": float(group.ppl_score.iloc[0]),
                     "models": len(group), "quality_pass_models": len(quality), "sequence_pass": sequence_pass,
                     "graphpep_median": median, "graphpep_mean": float(group.binding_score.mean()),
                     "graphpep_sd": float(group.binding_score.std(ddof=1)),
                     "graphpep_min": float(group.binding_score.min()), "graphpep_max": float(group.binding_score.max()),
                     "quality_pass_graphpep_median": float(quality.binding_score.median()),
                     "median_gain_vs_parent": parent_median - median, "matched_model_mean_gain": float(gain.mean()),
                     "models_improving_vs_parent": improvements, "predefined_improvement": improved,
                     "strict_three_model_improvement": bool(improved and improvements == 3),
                     "iptm_mean": float(group.iptm.mean()), "peptide_plddt_mean": float(group.peptide_plddt.mean()),
                     "interface_pae_median": float(group.interface_pae.median()),
                     "contact_jaccard_mean": float(group.contact_jaccard_other_models.mean()),
                     "native_contact_recall_mean": float(group.native_contact_recall_5A.mean()),
                     "peptide_CA_RMSD_mean_A": float(group.peptide_CA_RMSD_after_receptor_alignment_A.mean()),
                     "clash_models": int((group.clashing_atompairs_under_1p5A > 0).sum())})
    sequences = pd.DataFrame(rows)
    arm_rows = []
    for arm, group in sequences[sequences.arm != "parent"].groupby("arm"):
        models = frame[frame.arm == arm]
        arm_rows.append({"arm": arm, "n_children": len(group), "quality_pass_sequences": int(group.sequence_pass.sum()),
                         "quality_pass_models": int(models.structure_pass.sum()), "n_models": len(models),
                         "children_improved": int(group.predefined_improvement.sum()),
                         "children_improved_all_3_models": int(group.strict_three_model_improvement.sum()),
                         "mean_child_median_score": float(group.graphpep_median.mean()),
                         "mean_child_median_gain_vs_parent": float(group.median_gain_vs_parent.mean()),
                         "quality_pass_mean_child_median_score": float(group.loc[group.sequence_pass, "graphpep_median"].mean()),
                         "mean_contact_jaccard": float(group.contact_jaccard_mean.mean()),
                         "mean_native_contact_recall": float(group.native_contact_recall_mean.mean()),
                         "mean_peptide_CA_RMSD_A": float(group.peptide_CA_RMSD_mean_A.mean())})
    return sequences, pd.DataFrame(arm_rows)


def clean_json(value):
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_report(frame: pd.DataFrame, sequences: pd.DataFrame, arms: pd.DataFrame, config: dict) -> str:
    parent = sequences[sequences.arm == "parent"].iloc[0]
    guided = arms[arms.arm == "guided"].iloc[0]
    random = arms[arms.arm == "random"].iloc[0]
    if not parent.sequence_pass:
        conclusion = "本轮父代未通过预先规定的结构门槛，反馈策略的比较结论不确定。"
    elif guided.children_improved > random.children_improved:
        conclusion = "反馈组达到预设改善标准的子代更多，得到一个值得复验的计算信号；单父代、每组3条不足以证明普遍优势。"
    elif guided.children_improved == random.children_improved:
        conclusion = "两组达到预设改善标准的子代数量相同，本轮未显示反馈选位的明确优势。"
    else:
        conclusion = "随机组达到预设改善标准的子代更多，本轮不支持这条反馈选位规则优于随机选位。"
    lines = ["# MDM2–p53局部反馈闭环：一轮小规模对照实验", "", f"**结果：{conclusion}**", "",
             "## 实际完成", "", "- 1条实验结构来源父代，3条GraphPep反馈指导子代，3条随机mask对照。",
             "- 每条序列3个AFM v3参数模型，同一seed；共21个复合物预测及GraphPep评价。三个模型不是三个独立生物学重复。",
             "- 使用现有本地PepMLM、AFM和GraphPep。MSA复用既有受体缓存，本轮未请求公共MSA服务。未训练模型。", "",
             "## 预先固定的规则", "",
             "父代为1YCR中有坐标的MDM2 25–109和p53 17–29，肽序列为`ETFSDLWKLLPEN`。这与完整表达构建体不同，结果限于该截短构建体。",
             "GraphPep对实验父代的接触置信度按位点求均值，在有接触记录的位置中选最低两位：肽第4/8位，即p53 S20/K24。无接触的第5位不视为低分；随机组使用同一候选位置集合。",
             "两组每条子代都替换2个残基；PepMLM单次前向，标准20种氨基酸中排除父代残基后top-k=3、温度1，按概率采样。组内不重复；没有PPL或结构分数预筛选。未人工保护已知锚点。",
             "父代和子代使用完全相同的受体MSA，肽仅提供自身序列，无肽同源MSA、无模板。AFM模型1/2/3，seed 20261005，3 recycles且不提前停止，MSA 256/512，无relax。",
             "模型结构门槛：ipTM≥0.70，肽均值pLDDT≥70，界面PAE中位数≤10 Å，有跨链接触且最短重原子距离≥1.5 Å。序列至少2/3模型通过。这些是本轮探索性阈值。",
             "改善标准：父代及子代通过序列门槛，子代GraphPep中位数比预测父代降低>0.05，且至少2/3相同参数模型评分更好。另列出全部3模型改善。实验父代分数没有用作子代排名基线。", "",
             "## 每条序列结果", "", "GraphPep原始分数越负越好；该分数不等于亲和力。增益=预测父代中位数−子代中位数，正值为改善。", "",
             "| ID | 序列 | PPL | GraphPep中位数 | 对父代增益 | 结构通过 | 同模型改善 | 达到预设标准 |",
             "|---|---|---:|---:|---:|---:|---:|---|"]
    order = ["PARENT", "GUIDED_01", "GUIDED_02", "GUIDED_03", "RANDOM_01", "RANDOM_02", "RANDOM_03"]
    for ident in order:
        row = sequences.set_index("peptide_id").loc[ident]
        lines.append(f"| {ident} | `{row.sequence}` | {row.ppl_score:.2f} | {row.graphpep_median:.4f} | {row.median_gain_vs_parent:+.4f} | {row.quality_pass_models}/3 | {row.models_improving_vs_parent}/3 | {'是' if row.predefined_improvement else '否'} |")
    lines += ["", "## 两组比较", "", "| 组别 | 通过结构门槛的子代 | 达到改善标准 | 9个结构中通过数 | 子代中位数的平均增益 |",
              "|---|---:|---:|---:|---:|"]
    for _, row in arms.iterrows():
        lines.append(f"| {'反馈组' if row.arm == 'guided' else '随机组'} | {row.quality_pass_sequences}/3 | {row.children_improved}/3 | {row.quality_pass_models}/9 | {row.mean_child_median_gain_vs_parent:+.4f} |")
    lines += ["", "![模型评分与质量](analysis/model_scores.png)", "",
              "## 解释范围与下一步", "",
              "这次完成了实验结构反馈→选mask位置→PepMLM生成→重新预测复合物→GraphPep评价的一轮回路。反馈来源于实验父代，因此它验证的是有实验结构起点的局部优化，尚未证明预测结构反馈具有同样效果。",
              "GraphPep按超过0.5/0.7/0.9阈值的接触数量聚合，更多接触可能改善分数；低均值接触置信度不等于低能量贡献。需同时查看碰撞、PAE、接触稳定性及结构图。",
              "原生接触恢复率按残基位置计算，子代允许改变氨基酸；它与受体对齐后的肽CA RMSD都是姿态保留诊断，不能证明亲和力改善，也不排除其他合理结合模式。",
              "一条父代、每组3条序列是工程可行性试验。MDM2–p53与模型训练数据是否重叠尚未排除；结果不代表独立泛化，也不能推广为NCAM1有效。",
              "后续优先复验本轮候选的独立随机种子，再增加有实验结合数据的父代与突变体；要验证亲和力提高，需要结合实验。若扩展预测反馈闭环，应单独设置预测结构反馈选位组。", "",
              "## 数据", "",
              "- `protocol.json`：生成前固定的协议、种子、输入哈希及阈值。",
              "- `candidates.csv` / `sampling_log.json`：全部提案、PPL、mask位点、抽样概率。",
              "- `structures/`：21个预测、置信度与完整折叠日志。",
              "- `analysis/model_results.csv`：逐模型结构质量、GraphPep分数及反馈路径。",
              "- `analysis/sequence_summary.csv` / `analysis/arm_summary.csv`：逐序列、逐组比较。",
              "- `analysis/verification.json`：序列、MSA一致性、模型计数、评分公式和源文件哈希核验。",
              "- `analysis/graphpep/`：逐复合物输入、日志、接触与残基反馈。", "",
              "结构来源：[RCSB PDB 1YCR](https://www.rcsb.org/structure/1YCR)。"]
    (OUT / "report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return conclusion


def plot_scores(frame: pd.DataFrame, sequences: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = ["PARENT", "GUIDED_01", "GUIDED_02", "GUIDED_03", "RANDOM_01", "RANDOM_02", "RANDOM_03"]
    colors = {"parent": "#555555", "guided": "#24755c", "random": "#c77d31"}
    fig, ax = plt.subplots(figsize=(10, 5))
    for i, ident in enumerate(order):
        group = frame[frame.peptide_id == ident].sort_values("model_number")
        color = colors[group.arm.iloc[0]]
        for offset, (_, row) in zip([-.15, 0, .15], group.iterrows()):
            ax.scatter(i + offset, row.binding_score, color=color, marker="o" if row.structure_pass else "x", s=65)
        median = float(group.binding_score.median())
        ax.plot([i - .25, i + .25], [median, median], color=color, linewidth=2)
    baseline = float(sequences.loc[sequences.arm == "parent", "graphpep_median"].iloc[0])
    ax.axhline(baseline, linestyle="--", color="#888888", linewidth=1, label="Predicted parent median")
    ax.set_xticks(range(7), order, rotation=25, ha="right")
    ax.set_ylabel("GraphPep raw score (lower is better)")
    ax.set_title("MDM2-p53: one parent, six proposals, three AFM models each")
    ax.text(.01, .02, "Circle: structure gate passed; cross: failed; horizontal segment: median", transform=ax.transAxes, fontsize=9)
    ax.legend(loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(ANALYSIS / "model_scores.png", dpi=180)
    plt.close(fig)


def main(summary_only: bool = False) -> None:
    started = time.perf_counter()
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    config = json.loads((OUT / "protocol.json").read_text())
    for name, expected in {**config["source_hashes"], **config["msa_hashes"]}.items():
        if pilot.sha256(Path(name)) != expected:
            raise ValueError(f"Frozen input changed: {name}")
    if pilot.sha256(OUT / "candidates.csv") != config["candidate_file_sha256"]:
        raise ValueError("Candidate manifest changed after generation")
    candidates = pd.read_csv(OUT / "candidates.csv", keep_default_na=False)
    for row in candidates[candidates.arm != "parent"].to_dict("records"):
        pilot.validate_child(config["parent_sequence"], row["sequence"], list(map(int, row["mask_positions"].split(";"))))
    afm_config = json.loads((OUT / "structures/config.json").read_text())
    for name, expected in {"num_models": 3, "num_recycles": 3, "random_seed": 20261005,
                           "num_seeds": 1, "use_templates": False, "num_relax": 0,
                           "model_order": [1, 2, 3], "max_seq": 256, "max_extra_seq": 512,
                           "recycle_early_stop_tolerance": 0}.items():
        if afm_config[name] != expected:
            raise ValueError(f"AFM config mismatch: {name}")
    metrics = collect_metrics(config, candidates)
    metrics.to_csv(ANALYSIS / "structure_metrics.csv", index=False)
    if summary_only:
        frame = pd.read_csv(ANALYSIS / "model_results.csv")
    else:
        metrics.to_csv(ANALYSIS / "graphpep_manifest.csv", index=False)
        pairs = gp.build_inputs(metrics, str(ANALYSIS / "graphpep"))
        scores = []
        for index, pair in enumerate(pairs, 1):
            print(f"GraphPep {index}/21: {pair['row']['peptide_id']} {pair['row']['model_id']}", flush=True)
            result = gp.score_pair(pair, gp.GP, gp.ROOT / "pipeline/graphpep_patches/resmap_compat.py", .7)
            result.update(verify_graphpep_score(result))
            scores.append(result)
            pd.DataFrame(scores).to_csv(ANALYSIS / "model_results.csv", index=False)
        frame = pd.DataFrame(scores)
    if len(frame) != 21 or len(frame.peptide_id.unique()) != 7:
        raise ValueError("Incomplete scoring outputs")
    for row in frame.to_dict("records"):
        verify_graphpep_score(row)
    sequences, arms = summarize(frame, config)
    sequences.to_csv(ANALYSIS / "sequence_summary.csv", index=False)
    arms.to_csv(ANALYSIS / "arm_summary.csv", index=False)
    plot_scores(frame, sequences)
    gp_hashes = {path: pilot.sha256(Path(path)) for path in config["native_graphpep_source_hashes_before"]}
    if gp_hashes != config["native_graphpep_source_hashes_before"]:
        raise ValueError("Native GraphPep source hashes changed")
    verification = {"n_sequences": 7, "n_predictions": 21, "n_graphpep_evaluations": 21,
                    "n_finite_graphpep_scores": int(frame.binding_score.notna().sum()),
                    "sequence_and_confidence_dimensions_verified": True, "mutation_boundaries_verified": True,
                    "common_receptor_msa_hash": config["msa_records"][0]["common_receptor_rows_sha256"],
                    "frozen_inputs_unchanged": True, "score_formulas_verified_with_rounding_bounds": True,
                    "native_graphpep_source_hashes_unchanged": True, "native_graphpep_source_hashes": gp_hashes,
                    "afm_config_verified": True, "analysis_seconds": time.perf_counter() - started,
                    "analysis_script_sha256": pilot.sha256(Path(__file__)),
                    "score_table_sha256": pilot.sha256(ANALYSIS / "model_results.csv")}
    pilot.dump_json(ANALYSIS / "verification.json", verification)
    conclusion = write_report(frame, sequences, arms, config)
    summary = {"state": "completed", "conclusion_zh": conclusion,
               "parent": sequences[sequences.arm == "parent"].to_dict("records")[0],
               "arms": arms.to_dict("records"), "sequences": sequences.to_dict("records"), "verification": verification}
    pilot.dump_json(OUT / "summary.json", clean_json(summary))
    print(json.dumps(clean_json({"conclusion": conclusion, "arms": summary["arms"]}), ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    main(args.summary_only)
