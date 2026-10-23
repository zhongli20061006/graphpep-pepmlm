#!/usr/bin/env python3
"""Run a frozen one-parent, anchor-protected feedback-versus-random control."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import re
import time

os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
os.environ.setdefault("HF_HUB_OFFLINE", "1")

import numpy as np
import pandas as pd
import torch
from transformers import AutoModelForMaskedLM, AutoTokenizer

import feedback_panel as panel

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/feedback_control_1ycr_20261006"
SEED = 20261010
N_CHILDREN = 10
TOP_K = 4
PROTECTED = {3, 7, 10}
PARENT_ID = "P53WT_1YCR"

# Reuse the already validated 1YCR definition and helpers without changing the
# prior panel's output directory or seed.
PARENT = next(row.copy() for row in panel.PARENTS if row["peptide_id"] == PARENT_ID)
panel.OUT = OUT
panel.SEED = SEED
panel.PARENTS = [PARENT]
panel.COMPARISON_PROTOCOL = {
    **panel.COMPARISON_PROTOCOL,
    "afm_seed": SEED,
    "relaxation": {
        "method": "ColabFold Amber restrained minimization",
        "max_iterations": 200,
        "max_outer_iterations": 3,
        "tolerance": 2.39,
        "stiffness": 10.0,
        "use_gpu": False,
    },
}

PLAN = ROOT / "docs/superpowers/plans/2026-10-06-mdm2-1ycr-matched-feedback-control.md"
NATIVE_GP_FILES = ["GraphPep.sh", "bin/core.py", "bin/data.py"]
QUALITY_GATE = {
    "iptm_min": 0.70,
    "peptide_plddt_min": 70.0,
    "interface_pae_max": 10.0,
    "minimum_distance_min": 1.5,
    "minimum_contacts": 1,
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def load_metadata() -> dict:
    path = OUT / "parent_metadata.json"
    if not path.is_file():
        raise FileNotFoundError(f"Run prepare first; missing {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def prepare() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metadata_path = OUT / "parent_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("parents", [{}])[0].get("peptide_id") != PARENT_ID:
            raise ValueError("Existing output belongs to a different parent")
        print(f"Parent inputs already prepared: {metadata.get('stage')}", flush=True)
        return
    panel.make_parent_inputs()


def score_parent() -> None:
    metadata = load_metadata()
    if metadata.get("stage") not in {"native_controls_ready", "parent_structures_scored_feedback_frozen"}:
        raise ValueError(f"Unexpected parent stage: {metadata.get('stage')}")
    panel.parent_structure_scores()
    metadata = load_metadata()
    panel.verify_parent_metrics(metadata)
    metrics = pd.read_csv(OUT / "parent_model_metrics.csv")
    passing = metrics[metrics.structure_pass.astype(bool)]
    if len(passing) < 2:
        raise ValueError(f"Only {len(passing)}/3 parent models pass the frozen quality gate")
    print(f"Parent AFM models passing fixed gate: {len(passing)}/3", flush=True)


def observed_pool(profile: list[dict]) -> list[int]:
    pool = [int(row["peptide_position"]) for row in profile
            if int(row["count"]) > 0 and int(row["peptide_position"]) not in PROTECTED]
    if len(pool) < 2:
        raise ValueError(f"Need at least two observed non-anchor sites; got {pool}")
    return sorted(pool)


def current_source_hashes() -> dict[str, str]:
    parent_native = panel.PRIOR / "POSCTRL_MDM2_P53_native.pdb"
    cached_msa = panel.PRIOR / "results/POSCTRL_MDM2_P53.a3m"
    feedback = OUT / f"{PARENT_ID}_predicted_ensemble_feedback.csv"
    predicted_scores = OUT / "parent_model_scores.csv"
    predicted_metrics = OUT / "parent_model_metrics.csv"
    sources = [Path(__file__), Path(panel.__file__), ROOT / "pipeline/shared_receptor_msa.py",
               ROOT / "pipeline/step1_generate.py", ROOT / "pipeline/relax_feedback_panel.py",
               PLAN, parent_native, cached_msa, feedback, predicted_scores, predicted_metrics]
    sources += [panel.gp.GP / name for name in NATIVE_GP_FILES]
    for name in ("config.json", "tokenizer_config.json", "special_tokens_map.json"):
        candidate = panel.MODEL / name
        if candidate.is_file():
            sources.append(candidate)
    unique = list(dict.fromkeys(sources))
    missing = [str(path) for path in unique if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing frozen protocol source(s): " + ", ".join(missing))
    return {str(path): digest(path) for path in unique}


def freeze_and_generate() -> None:
    metadata = load_metadata()
    if metadata.get("stage") != "parent_structures_scored_feedback_frozen":
        raise ValueError("Score and freeze parent AFM feedback before generating offspring")
    if not metadata["parents"] or metadata["parents"][0]["peptide_id"] != PARENT_ID:
        raise ValueError("Frozen parent metadata does not describe 1YCR")
    candidate_path = OUT / "candidates.csv"
    if candidate_path.exists():
        if metadata.get("candidate_sha256") == digest(candidate_path):
            print("Candidate set is already frozen; reusing it", flush=True)
            return
        raise ValueError("Candidate file exists but does not match parent metadata")

    parent = metadata["parents"][0]
    parent_sequence = parent["sequence"]
    if parent_sequence != "ETFSDLWKLLPEN":
        raise ValueError(f"Unexpected 1YCR peptide sequence: {parent_sequence}")
    parent_metrics = pd.read_csv(OUT / "parent_model_metrics.csv")
    passing_ids = set(parent_metrics.loc[parent_metrics.structure_pass.astype(bool), "model_id"].astype(str))
    if len(passing_ids) < 2:
        raise ValueError("At least two parent structures must pass before site selection")

    feedback_path = Path(parent["predicted_ensemble_profile_path"])
    profile = panel.profile(str(feedback_path), len(parent_sequence))
    pool = observed_pool(profile)
    ranked = sorted((row for row in profile if row["peptide_position"] in pool),
                    key=lambda row: (float(row["mean"]), int(row["peptide_position"])))
    guided_sites = sorted(int(row["peptide_position"]) for row in ranked[:2])
    if len(guided_sites) != 2 or set(guided_sites) & PROTECTED:
        raise ValueError(f"Invalid feedback-selected sites: {guided_sites}")

    possible_pairs = [tuple(pair) for pair in itertools.combinations(pool, 2)
                      if tuple(pair) != tuple(guided_sites)]
    if len(possible_pairs) < N_CHILDREN:
        raise ValueError("Not enough non-guided random site pairs for ten controls")
    position_rng = np.random.default_rng(SEED)
    chosen_indices = position_rng.choice(len(possible_pairs), size=N_CHILDREN, replace=False)
    random_pairs = [list(possible_pairs[int(index)]) for index in chosen_indices]

    source_hashes = current_source_hashes()
    if digest(PLAN) != source_hashes[str(PLAN)]:
        raise ValueError("Plan hash could not be frozen")
    protocol = {
        "state": "frozen_before_offspring_scoring",
        "parent_id": PARENT_ID,
        "parent_sequence": parent_sequence,
        "parent_pdb": "1YCR",
        "protected_sites": [{"position": p, "residue": parent_sequence[p - 1],
                              "p53_pdb_residue": p + 16}
                             for p in sorted(PROTECTED)],
        "feedback_source": "contact-count-weighted GraphPep residue feedback pooled over quality-passing AFM-v3 parent models",
        "feedback_model_ids": sorted(passing_ids),
        "feedback_profile_path": str(feedback_path),
        "feedback_profile_sha256": digest(feedback_path),
        "eligible_observed_nonanchor_positions": pool,
        "guided_positions": guided_sites,
        "random_position_pairs": random_pairs,
        "random_pair_sampling": "uniform without replacement over distinct 2-site pairs from the shared eligible pool; guided pair excluded",
        "children_per_arm": N_CHILDREN,
        "mutations_per_child": 2,
        "top_k": TOP_K,
        "top_k_rationale": "top-k 3 gives at most 9 unique joint sequences at the fixed guided pair; top-k 4 is the smallest shared value that permits 10 unique children",
        "temperature": 1.0,
        "ppl_filter": None,
        "generation_seed": SEED,
        "pepmlm_draw_seeds": {
            "predicted_guided": [SEED + 1000 + i for i in range(1, N_CHILDREN + 1)],
            "predicted_random": [SEED + 2000 + i for i in range(1, N_CHILDREN + 1)],
        },
        "afm": {"model_type": "alphafold2_multimer_v3", "models": [1, 2, 3],
                "random_seed": SEED, "recycles": 3, "recycle_early_stop_tolerance": 0,
                "max_seq": 256, "max_extra_seq": 512, "templates": False,
                "num_relax": 0, "unified_memory": False},
        "amber_relaxation": panel.COMPARISON_PROTOCOL["relaxation"],
        "quality_gate": QUALITY_GATE,
        "improvement_rule": {"minimum_common_passing_models": 2, "minimum_graphpep_gain": 0.05,
                              "minimum_matched_models_improved": 2,
                              "gain_definition": "parent GraphPep score minus child GraphPep score; positive is lower child score"},
        "source_hashes": source_hashes,
        "graphpep_source_hashes": {str(panel.gp.GP / name): digest(panel.gp.GP / name) for name in NATIVE_GP_FILES},
        "prior_parent_graphpep_metadata_sha256": digest(OUT / "parent_model_scores.csv"),
        "software_plan": str(PLAN),
        "source_snapshot": {"pepmlm_model_dir": str(panel.MODEL),
                            "pepmlm_config_sha256": source_hashes.get(str(panel.MODEL / "config.json"))},
    }
    write_json(OUT / "protocol.json", protocol)

    if not torch.cuda.is_available():
        raise RuntimeError("The validated local PepMLM workflow requires CUDA")
    device = "cuda"
    torch.manual_seed(SEED)
    tokenizer = AutoTokenizer.from_pretrained(panel.MODEL, local_files_only=True)
    model = AutoModelForMaskedLM.from_pretrained(
        panel.MODEL, torch_dtype=torch.float16, local_files_only=True).to(device).eval()
    canonical = [tokenizer.convert_tokens_to_ids(aa) for aa in panel.AMINO_ACIDS]
    if len(set(canonical)) != 20 or tokenizer.unk_token_id in canonical:
        raise ValueError("PepMLM canonical amino-acid token mapping is invalid")
    print(f"PepMLM loaded on {torch.cuda.get_device_name(0)}; guided sites={guided_sites}; eligible={pool}", flush=True)

    candidate_rows = [{"peptide_id": PARENT_ID, "parent_id": PARENT_ID, "arm": "parent",
                       "replicate": 0, "sequence": parent_sequence, "length": len(parent_sequence),
                       "mask_positions": "", "p53_positions": "", "mutations": "",
                       "mutation_count": 0, "ppl_score": math.nan, "sampling_seed": SEED}]
    proposals = []
    seen = {parent_sequence}
    arms = [("predicted_guided", [guided_sites] * N_CHILDREN,
             [SEED + 1000 + i for i in range(1, N_CHILDREN + 1)]),
            ("predicted_random", random_pairs,
             [SEED + 2000 + i for i in range(1, N_CHILDREN + 1)])]
    for arm, site_pairs, draw_seeds in arms:
        for replicate, (positions, draw_seed) in enumerate(zip(site_pairs, draw_seeds), 1):
            ident = f"{arm.upper()}_{PARENT_ID}_{replicate:02d}"
            child, trace = panel.draw_child(parent_sequence, positions, model, tokenizer,
                                            device, draw_seed, seen, top_k=TOP_K)
            mutations = ";".join(f"{parent_sequence[pos - 1]}{pos}{child[pos - 1]}" for pos in positions)
            candidate_rows.append({"peptide_id": ident, "parent_id": PARENT_ID, "arm": arm,
                                   "replicate": replicate, "sequence": child, "length": len(child),
                                   "mask_positions": ";".join(map(str, positions)),
                                   "p53_positions": ";".join(str(pos + 16) for pos in positions),
                                   "mutations": mutations, "mutation_count": len(positions),
                                   "ppl_score": math.nan, "sampling_seed": draw_seed})
            proposals.append({"peptide_id": ident, "arm": arm, "replicate": replicate,
                              "site_pool": pool, "selected_positions": positions,
                              "sampling_seed": draw_seed, "mutation_summary": mutations,
                              "sampling_trace": trace, "sequence": child})
            print(f"{ident}: {child} ({mutations})", flush=True)

    frame = pd.DataFrame(candidate_rows)
    target_ids = torch.tensor(tokenizer(panel.TARGET)["input_ids"][:-1])
    ppl_values = panel.step1_generate.pseudo_ppl_chunk(
        model, tokenizer, target_ids, frame.sequence.tolist(), device, len(parent_sequence), 16)
    frame["ppl_score"] = ppl_values
    frame.to_csv(candidate_path, index=False)
    (OUT / "candidate_complexes.fasta").write_text(
        "".join(f">{row.peptide_id}\n{panel.TARGET}:{row.sequence}\n"
                for row in frame.itertuples()), encoding="utf-8")
    write_json(OUT / "sampling_log.json", {"seed": SEED, "top_k": TOP_K,
                                             "temperature": 1.0, "proposals": proposals})

    cached = (panel.PRIOR / "results/POSCTRL_MDM2_P53.a3m").read_text(encoding="utf-8")
    msa_dir = OUT / "candidate_msa"
    offspring_dir = OUT / "offspring_msa"
    msa_dir.mkdir(parents=True, exist_ok=True)
    offspring_dir.mkdir(parents=True, exist_ok=True)
    msa_hashes, msa_records = {}, []
    for row in candidate_rows:
        msa, stats = panel.make_shared_msa(cached, panel.TARGET, parent_sequence, row["sequence"])
        path = msa_dir / f"{row['peptide_id']}.a3m"
        path.write_text(msa, encoding="utf-8")
        msa_hashes[str(path)] = digest(path)
        msa_records.append({"peptide_id": row["peptide_id"], **stats})
        if row["arm"] != "parent":
            (offspring_dir / path.name).write_text(msa, encoding="utf-8")
    receptor_hashes = {row["receptor_alignment_sha256"] for row in msa_records}
    if len(receptor_hashes) != 1:
        raise ValueError("Receptor MSA rows differ across parent and offspring")
    frame.to_csv(OUT / "candidates.csv", index=False)
    metadata.update(stage="controlled_offspring_and_msas_ready",
                    candidate_sha256=digest(candidate_path), protocol_sha256=digest(OUT / "protocol.json"),
                    candidate_msa_hashes=msa_hashes, msa_records=msa_records,
                    receptor_alignment_sha256=next(iter(receptor_hashes)),
                    pepmlm_device=torch.cuda.get_device_name(0),
                    software={"torch": torch.__version__, "numpy": np.__version__})
    write_json(OUT / "parent_metadata.json", metadata)
    print(f"Frozen {len(frame) - 1} unique offspring in two arms", flush=True)


def candidate_metrics(metadata: dict, candidates: pd.DataFrame) -> pd.DataFrame:
    parent = metadata["parents"][0]
    parent_spec = {**parent, "gate": QUALITY_GATE}
    native = Path(parent["native_path"])
    native_contacts = panel.collect_native_contact_positions(parent)
    rows = []
    for candidate in candidates.to_dict("records"):
        folder = OUT / "parent_relaxed" if candidate["arm"] == "parent" else OUT / "relaxed_structures"
        paths = panel.predictions(folder, candidate)
        for path in paths:
            rows.append(panel.metrics(parent_spec, candidate, path, native, native_contacts))
    return pd.DataFrame(rows)


def analyze() -> None:
    metadata = load_metadata()
    if metadata.get("stage") != "controlled_offspring_and_msas_ready":
        raise ValueError("Generate and freeze controlled candidates before scoring")
    if digest(OUT / "candidates.csv") != metadata["candidate_sha256"]:
        raise ValueError("Frozen candidates changed")
    if digest(OUT / "protocol.json") != metadata["protocol_sha256"]:
        raise ValueError("Frozen experiment protocol changed")
    protocol = json.loads((OUT / "protocol.json").read_text(encoding="utf-8"))
    for raw, expected in protocol["source_hashes"].items():
        if digest(Path(raw)) != expected:
            raise ValueError(f"Frozen source changed: {raw}")
    for raw, expected in metadata["candidate_msa_hashes"].items():
        if digest(Path(raw)) != expected:
            raise ValueError(f"Frozen MSA changed: {raw}")

    candidates = pd.read_csv(OUT / "candidates.csv", keep_default_na=False)
    children = candidates[candidates.arm != "parent"].copy()
    if len(children) != 2 * N_CHILDREN or children.sequence.nunique() != len(children):
        raise ValueError("Expected twenty unique offspring")
    if set(children.arm.value_counts().to_dict().items()) != {
            ("predicted_guided", N_CHILDREN), ("predicted_random", N_CHILDREN)}:
        raise ValueError("Each experiment arm must have exactly ten children")
    for row in children.to_dict("records"):
        positions = [int(value) for value in row["mask_positions"].split(";")]
        if len(positions) != 2 or len(set(positions)) != 2 or set(positions) & PROTECTED:
            raise ValueError(f"Invalid or anchor-mutating mask: {row['peptide_id']} {positions}")
        panel.common.validate_child(protocol["parent_sequence"], row["sequence"], positions)

    for raw, expected in protocol["graphpep_source_hashes"].items():
        if digest(Path(raw)) != expected:
            raise ValueError(f"GraphPep source changed: {raw}")

    analysis = OUT / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    model_metrics = candidate_metrics(metadata, candidates)
    model_metrics.to_csv(analysis / "model_metrics.csv", index=False)
    offspring_manifest = model_metrics[model_metrics.arm != "parent"].rename(
        columns={"model_path": "complex_path"})
    offspring_scores = panel.run_graphpep(offspring_manifest, analysis / "offspring_graphpep")
    parent_scores = pd.read_csv(OUT / "parent_model_scores.csv")
    score_rows = []
    for metric in model_metrics.to_dict("records"):
        source = parent_scores if metric["arm"] == "parent" else offspring_scores
        match = source[source.complex_path.astype(str) == str(metric["model_path"])]
        if len(match) != 1:
            raise ValueError(f"Expected one GraphPep record for {metric['model_path']}")
        row = {**metric, **match.iloc[0].to_dict()}
        panel.score_bounds(row)
        score_rows.append(row)
    results = pd.DataFrame(score_rows)
    if len(results) != 63:
        raise ValueError(f"Expected 63 predicted model records; found {len(results)}")
    results.to_csv(analysis / "model_results.csv", index=False)

    parent_rows = results[results.arm == "parent"].set_index("model_number")
    sequence_rows = []
    for peptide_id, group in results[results.arm != "parent"].groupby("peptide_id", sort=False):
        group = group.set_index("model_number").sort_index()
        paired = group.structure_pass.astype(bool) & parent_rows.loc[group.index].structure_pass.astype(bool)
        child_pass = group[group.structure_pass.astype(bool)]
        qualified_child = group[paired]
        qualified_parent = parent_rows.loc[qualified_child.index]
        gains = qualified_parent.binding_score.astype(float) - qualified_child.binding_score.astype(float)
        median_gain = float(gains.median()) if len(gains) >= 2 else math.nan
        improved_models = int((gains > 0).sum())
        passes_rule = bool(len(gains) >= 2 and median_gain > 0.05 and improved_models >= 2)
        candidate = children[children.peptide_id == peptide_id].iloc[0]
        sequence_rows.append({
            "peptide_id": peptide_id, "arm": candidate.arm, "replicate": int(candidate.replicate),
            "sequence": candidate.sequence, "mask_positions": candidate.mask_positions,
            "mutations": candidate.mutations, "ppl_score_diagnostic_only": float(candidate.ppl_score),
            "quality_pass_models": int(child_pass.shape[0]),
            "common_parent_child_passing_models": int(len(gains)),
            "graphpep_median_all_models": float(group.binding_score.median()),
            "parent_graphpep_median": float(parent_rows.binding_score.median()),
            "median_matched_gain": median_gain, "matched_models_improved": improved_models,
            "predefined_improvement": passes_rule,
        })
    sequence_summary = pd.DataFrame(sequence_rows)
    sequence_summary.to_csv(analysis / "sequence_summary.csv", index=False)

    arm_rows = []
    for arm, group in sequence_summary.groupby("arm", sort=False):
        finite = group.median_matched_gain[np.isfinite(group.median_matched_gain)]
        model_group = results[results.arm == arm]
        arm_rows.append({"arm": arm, "n_children": len(group),
                         "children_with_comparable_scores": int(np.isfinite(group.median_matched_gain).sum()),
                         "structure_pass_models": int(model_group.structure_pass.sum()),
                         "no_contact_models": int((model_group.graphpep_status == "no_contact").sum()),
                         "children_improved": int(group.predefined_improvement.sum()),
                         "mean_child_median_graphpep_gain": float(finite.mean()) if len(finite) else math.nan,
                         "median_child_median_graphpep_gain": float(finite.median()) if len(finite) else math.nan})
    arm_summary = pd.DataFrame(arm_rows)
    arm_summary.to_csv(analysis / "arm_summary.csv", index=False)

    records = candidates[candidates.arm != "parent"].merge(
        sequence_summary, on=["peptide_id", "arm", "replicate", "sequence", "mask_positions", "mutations"],
        suffixes=("", "_summary"), validate="one_to_one")
    feedback_paths = results[results.arm != "parent"].groupby("peptide_id").contact_feedback_path.apply(
        lambda paths: ";".join(paths.astype(str))).rename("residue_feedback_paths")
    structure_paths = results[results.arm != "parent"].groupby("peptide_id").model_path.apply(
        lambda paths: ";".join(paths.astype(str))).rename("structure_paths")
    records = records.merge(feedback_paths, on="peptide_id").merge(structure_paths, on="peptide_id")
    records["parent_source_pdb"] = "1YCR"
    records["structure_prediction_method"] = "ColabFold AFM-v3 models 1-3 + restrained Amber minimization"
    records.to_csv(OUT / "experiment_record.csv", index=False)

    parent_model_metrics = pd.read_csv(OUT / "parent_model_metrics.csv")
    passing_parents = int(parent_model_metrics.structure_pass.sum())
    receptor_hashes = {row["receptor_alignment_sha256"] for row in metadata["msa_records"]}
    score_status_ok = results.graphpep_status.isin(["ok", "no_contact"]).all()
    if passing_parents < 2 or len(receptor_hashes) != 1 or not score_status_ok:
        raise ValueError("Parent, MSA, or GraphPep integrity gate failed")
    gp_after = {raw: digest(Path(raw)) for raw in protocol["graphpep_source_hashes"]}
    if gp_after != protocol["graphpep_source_hashes"]:
        raise ValueError("GraphPep source hashes changed during scoring")
    verification = {
        "experiment": "MDM2 1YCR matched predicted-feedback and random-site control",
        "parent_count": 1, "parent_models": 3, "parent_models_passing": passing_parents,
        "offspring_count": len(children), "offspring_per_arm": N_CHILDREN,
        "offspring_structure_models": int((results.arm != "parent").sum()),
        "all_predicted_model_records": len(results),
        "graphpep_ok_records": int((results.graphpep_status == "ok").sum()),
        "graphpep_no_contact_records": int((results.graphpep_status == "no_contact").sum()),
        "finite_graphpep_scores": int(results.binding_score.notna().sum()),
        "all_sequences_unique": children.sequence.nunique() == len(children),
        "all_children_have_exactly_two_mutations": bool((children.mutation_count.astype(int) == 2).all()),
        "no_protected_anchor_mutations": bool(all(
            not (set(int(x) for x in row.mask_positions.split(";")) & PROTECTED)
            for row in children.itertuples())),
        "common_receptor_msa_rows": len(receptor_hashes) == 1,
        "graphpep_source_hashes_unchanged": gp_after == protocol["graphpep_source_hashes"],
        "protocol_hash_unchanged": digest(OUT / "protocol.json") == metadata["protocol_sha256"],
        "seed": SEED,
    }
    write_json(analysis / "verification.json", verification)
    write_report(protocol, arm_summary, verification)
    print(arm_summary.to_string(index=False), flush=True)
    print(f"Wrote {len(records)} peptide records and 63 model records", flush=True)


def write_report(protocol: dict, arm_summary: pd.DataFrame, verification: dict) -> None:
    def value(arm: str, column: str):
        found = arm_summary[arm_summary.arm == arm]
        return found.iloc[0][column] if len(found) else math.nan

    lines = [
        "# 1YCR 父代：预测反馈选位与匹配随机对照",
        "",
        "本轮只评估已解析完整肽序列的 1YCR MDM2–p53 复合物。两组各 10 条子代，每条替换 2 个非锚点位点；F3、W7、L10 固定不变。预测反馈来自父代 AFM 参数模型中通过结构质量门槛模型的 GraphPep 接触置信度。",
        "",
        f"父代预测反馈选位：`{protocol['guided_positions']}`；共同可选位点池：`{protocol['eligible_observed_nonanchor_positions']}`。随机组选位对从同一位点池均匀抽取，且不重复抽到引导组的位点对。",
        "",
        f"为保证固定引导位点对仍可生成 10 条不同序列，两组均使用 top-k={protocol['top_k']}。top-k=3 在两个固定 mask 位点上至多产生 9 种序列；因此 top-k=4 是满足 10 条 unique 子代的最小共同设置。PPL 仅作记录，没有用于筛选。",
        "",
        "结构预测使用同一个 AFM-v3 seed 和参数模型 1/2/3；父代和子代都经过相同的 Amber 有约束几何松弛。GraphPep 按模型编号与父代配对。改善定义沿用预先冻结的门槛：至少两个父子配对模型通过结构筛选，GraphPep 配对中位增益大于 0.05，且至少两个模型的子代分数更低。",
        "",
        "## 组间结果",
        "",
        "|选位组|子代数|有可比评分|结构通过模型|无接触模型|达标子代|平均子代中位增益|中位子代中位增益|",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for arm, label in [("predicted_guided", "AFM反馈引导"), ("predicted_random", "匹配随机")]:
        lines.append(f"|{label}|{int(value(arm, 'n_children'))}|{int(value(arm, 'children_with_comparable_scores'))}|"
                     f"{int(value(arm, 'structure_pass_models'))}|{int(value(arm, 'no_contact_models'))}|"
                     f"{int(value(arm, 'children_improved'))}/10|{value(arm, 'mean_child_median_graphpep_gain'):.4f}|"
                     f"{value(arm, 'median_child_median_graphpep_gain'):.4f}|")
    lines += [
        "",
        "## 解释范围",
        "",
        "这是一组单父代、小规模的计算对照。AFM 参数模型不是生物学重复，十条子代也不足以证明普遍规律；结果只用于判断在保护已知锚点并匹配位点池后，预测反馈选位是否显示更一致的 GraphPep 代理评分变化。GraphPep 分数不是 kcal/mol 的 ΔG(binding)，也不等于实测 Kd。",
        "",
        f"完整性记录：父代结构通过 {verification['parent_models_passing']}/3；子代 {verification['offspring_count']} 条、AFM 子代结构模型 {verification['offspring_structure_models']} 个；GraphPep 有效接触评分 {verification['graphpep_ok_records']} 个、无接触记录 {verification['graphpep_no_contact_records']} 个。",
        "",
        "逐模型数据见 `analysis/model_results.csv`，逐肽和组级数据见 `analysis/sequence_summary.csv`、`analysis/arm_summary.csv`，冻结设置及哈希见 `protocol.json` 和 `analysis/verification.json`。",
        "",
    ]
    (OUT / "report_zh.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "score-parent", "generate", "analyze"])
    stage = parser.parse_args().stage
    started = time.perf_counter()
    if stage == "prepare":
        prepare()
    elif stage == "score-parent":
        score_parent()
    elif stage == "generate":
        freeze_and_generate()
    else:
        analyze()
    print(f"Stage {stage} completed in {time.perf_counter() - started:.1f}s", flush=True)


if __name__ == "__main__":
    main()
