#!/usr/bin/env python3
"""Generate a frozen six-proposal local-regeneration pilot around resolved 1YCR."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / "output/experiment0/expanded_validation"
OUT = ROOT / "output/feedback_pilot_mdm2_20261005"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
MODEL = Path("/home/g/pepmlm_pipeline/models/PepMLM-650M")
GP = Path("/home/g/GraphPep_v1.3")
SEED = 20261005


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def select_positions(rows: list[dict], number: int) -> tuple[list[int], list[int]]:
    """Missing contact information is ineligible, rather than a zero score."""
    observed = []
    for row in rows:
        if int(row["count"]) <= 0:
            continue
        value = float(row["mean"])
        if math.isfinite(value):
            if not 0 <= value <= 1:
                raise ValueError("GraphPep confidence outside [0,1]")
            observed.append((int(row["peptide_position"]), value))
    if len(observed) < number or len({p for p, _ in observed}) != len(observed):
        raise ValueError("Insufficient or duplicate observed positions")
    return sorted(p for p, _ in observed), [p for p, _ in sorted(observed, key=lambda x: (x[1], x[0]))[:number]]


def parse_a3m(text: str) -> list[tuple[str, str]]:
    records = []
    header, seq = None, []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq)))
            header, seq = line[1:], []
        else:
            if header is None:
                raise ValueError("Sequence before A3M header")
            seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def aligned(sequence: str) -> str:
    return "".join(c for c in sequence if not c.islower())


def make_shared_msa(text: str, receptor: str, parent: str, child: str) -> tuple[str, dict]:
    """Retain receptor-only rows, preserving insertions and dropping all peptide hits."""
    nr, np = len(receptor), len(parent)
    records = parse_a3m(text)
    if len(child) != np or not records or aligned(records[0][1]) != receptor + parent:
        raise ValueError("Cached MSA query does not match the fixed parent construct")
    common = []
    for header, sequence in records[1:]:
        clean = aligned(sequence)
        if len(clean) != nr + np:
            raise ValueError("Unexpected A3M aligned length")
        if any(c != "-" for c in clean[:nr]) and set(clean[nr:]) == {"-"}:
            common.append((header, sequence))
    if not common:
        raise ValueError("No receptor-only MSA rows")
    body = "".join(f">{h}\n{s}\n" for h, s in common)
    result = f"#{nr},{np}\t1,1\n>101\t102\n{receptor}{child}\n" + body + f">102\n{'-' * nr}{child}\n"
    stats = {"source_records": len(records), "receptor_only_records": len(common),
             "discarded_non_receptor_records": len(records) - 1 - len(common),
             "common_receptor_rows_sha256": hashlib.sha256(body.encode()).hexdigest()}
    return result, stats


def validate_child(parent: str, child: str, positions: list[int]) -> None:
    if len(child) != len(parent) or set(child) - set(AMINO_ACIDS):
        raise ValueError("Invalid child sequence or length")
    observed = {i + 1 for i, (a, b) in enumerate(zip(parent, child)) if a != b}
    if observed != set(positions) or len(positions) != len(set(positions)):
        raise ValueError(f"Mutation boundary mismatch: {observed} != {positions}")


def generate() -> None:
    import numpy as np
    import pandas as pd
    import torch
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    from step1_generate import pseudo_ppl_chunk

    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "candidates.csv").exists():
        raise FileExistsError("Candidates already frozen; reuse them instead of regenerating")
    metadata = json.loads((PRIOR / "positive_control.json").read_text())
    receptor, parent = [metadata["sequences_from_resolved_coordinates"][c] for c in ("A", "B")]
    native = PRIOR / "POSCTRL_MDM2_P53_native.pdb"
    if sha256(native) != metadata["native_complex_sha256"]:
        raise ValueError("Native parent hash changed")
    feedback = PRIOR / "analysis/graphpep/graphpep_inputs/0015_POSCTRL_MDM2_P53_1YCR_native/residue_feedback.csv"
    rows = list(csv.DictReader(feedback.open()))
    eligible, guided = select_positions(rows, 2)
    position_rng = np.random.default_rng(SEED)
    proposals = []
    for replicate in range(1, 4):
        random_positions = sorted(position_rng.choice(eligible, size=2, replace=False).tolist())
        for arm, positions in [("guided", sorted(guided)), ("random", random_positions)]:
            proposals.append({"peptide_id": f"{arm.upper()}_{replicate:02d}", "arm": arm,
                              "replicate": replicate, "mask_positions": positions, "seed": SEED + replicate})
    source_hashes = {str(path): sha256(path) for path in [native, feedback, PRIOR / "results/POSCTRL_MDM2_P53.a3m",
                                                       Path(__file__), ROOT / "pipeline/step1_generate.py"]}
    native_gp_hashes = {str(GP / name): sha256(GP / name) for name in ["GraphPep.sh", "bin/core.py", "bin/data.py"]}
    config = {"state": "protocol_frozen", "parent_sequence": parent, "receptor_sequence": receptor,
              "source_pdb": "1YCR", "p53_start": 17, "mdm2_start": 25,
              "feedback_source": str(feedback), "eligible_positions": eligible, "guided_positions": guided,
              "proposals": proposals, "top_k": 3, "temperature": 1., "substitutions": 2,
              "parent_residue_excluded_before_top_k": True, "model_dir": str(MODEL), "dtype": "fp16",
              "within_arm_duplicates": "joint alternatives renormalized after excluding already sampled sequences",
              "source_hashes": source_hashes, "native_graphpep_source_hashes_before": native_gp_hashes,
              "afm_models": [1, 2, 3], "afm_seed": SEED, "afm_recycles": 3,
              "afm_recycle_early_stop_tolerance": 0, "max_seq": 256, "max_extra_seq": 512,
              "templates": False, "relaxation": False, "unified_memory": False,
              "model_gate": {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10.,
                             "minimum_interchain_heavyatom_distance_min": 1.5, "minimum_contacts": 1},
              "sequence_gate": "at least 2 of 3 models pass", "minimum_score_gain": .05,
              "protocol_file": str(ROOT / "docs/superpowers/plans/2026-10-05-mdm2-feedback-pilot.md")}
    dump_json(OUT / "protocol.json", config)
    print(f"Frozen guided positions: {guided}; random positions: {[p['mask_positions'] for p in proposals if p['arm']=='random']}", flush=True)
    if not torch.cuda.is_available():
        raise RuntimeError("Expected the validated local CUDA device")
    torch.manual_seed(SEED)
    tok = AutoTokenizer.from_pretrained(MODEL, local_files_only=True)
    model = AutoModelForMaskedLM.from_pretrained(MODEL, torch_dtype=torch.float16, local_files_only=True).cuda().eval()
    print(f"PepMLM loaded on {torch.cuda.get_device_name(0)}", flush=True)
    canonical_ids = [tok.convert_tokens_to_ids(a) for a in AMINO_ACIDS]
    if len(set(canonical_ids)) != 20 or tok.unk_token_id in canonical_ids:
        raise ValueError("Canonical amino acid token mapping invalid")
    seen = {"guided": set(), "random": set()}
    candidates = [{"peptide_id": "PARENT", "arm": "parent", "replicate": 0, "sequence": parent,
                   "length": len(parent), "mask_positions": "", "p53_positions": "", "mutations": "",
                   "mutation_count": 0, "sampling_seed": SEED}]
    sampling_log = []
    for proposal in proposals:
        positions = proposal["mask_positions"]
        masked_peptide = "".join("<mask>" if i in positions else a for i, a in enumerate(parent, 1))
        encoded = tok(receptor + masked_peptide, return_tensors="pt").to("cuda")
        actual = (encoded.input_ids[0] == tok.mask_token_id).nonzero(as_tuple=True)[0].tolist()
        expected = [len(receptor) + p for p in positions]
        if actual != expected:
            raise ValueError(f"Token offsets mismatch: {actual} != {expected}")
        with torch.inference_mode():
            logits = model(**encoded).logits[0, actual].float().cpu()
        alternatives = []
        for i, position in enumerate(positions):
            values = logits[i, canonical_ids].clone()
            values[AMINO_ACIDS.index(parent[position - 1])] = -torch.inf
            values, indices = values.topk(3)
            probs = torch.softmax(values, dim=0).tolist()
            alternatives.append([(AMINO_ACIDS[j], p) for j, p in zip(indices.tolist(), probs)])
        joint = []
        for choice in itertools.product(*alternatives):
            child = list(parent)
            for pos, (aa, _) in zip(positions, choice):
                child[pos - 1] = aa
            sequence = "".join(child)
            if sequence not in seen[proposal["arm"]]:
                joint.append((sequence, math.prod(p for _, p in choice)))
        if not joint:
            raise ValueError("No unseen child alternatives remain")
        probs = np.asarray([p for _, p in joint], dtype=float)
        probs /= probs.sum()
        draw_rng = np.random.default_rng(proposal["seed"])
        selection = int(draw_rng.choice(len(joint), p=probs))
        child = joint[selection][0]
        validate_child(parent, child, positions)
        seen[proposal["arm"]].add(child)
        record = {"peptide_id": proposal["peptide_id"], "arm": proposal["arm"],
                  "replicate": proposal["replicate"], "sequence": child, "length": len(child),
                  "mask_positions": ";".join(map(str, positions)),
                  "p53_positions": ";".join(str(p + 16) for p in positions),
                  "mutations": ";".join(f"{parent[p-1]}{p}{child[p-1]}" for p in positions),
                  "mutation_count": len(positions), "sampling_seed": proposal["seed"]}
        candidates.append(record)
        sampling_log.append({**proposal, "alternatives": alternatives, "joint_available": joint,
                             "conditional_probabilities": probs.tolist(), "selected_sequence": child,
                             "selected_probability": float(probs[selection])})
        print(f"{record['peptide_id']}: {child} ({record['mutations']})", flush=True)
    target_ids = torch.tensor(tok(receptor)["input_ids"][:-1])
    ppls = pseudo_ppl_chunk(model, tok, target_ids, [r["sequence"] for r in candidates], "cuda", len(parent), 7)
    for record, ppl in zip(candidates, ppls):
        record["ppl_score"] = float(ppl)
    pd.DataFrame(candidates).to_csv(OUT / "candidates.csv", index=False)
    dump_json(OUT / "sampling_log.json", {"proposals": sampling_log})
    msa_dir = OUT / "msa_inputs"
    msa_dir.mkdir(exist_ok=True)
    cached = (PRIOR / "results/POSCTRL_MDM2_P53.a3m").read_text()
    msa_stats, msa_hashes = [], {}
    for record in candidates:
        msa, stats = make_shared_msa(cached, receptor, parent, record["sequence"])
        path = msa_dir / f"{record['peptide_id']}.a3m"
        path.write_text(msa)
        msa_hashes[str(path)] = sha256(path)
        msa_stats.append({"peptide_id": record["peptide_id"], **stats})
    if len({r["common_receptor_rows_sha256"] for r in msa_stats}) != 1:
        raise ValueError("MSA receptor input differs between arms")
    (OUT / "complexes.fasta").write_text("".join(f">{r['peptide_id']}\n{receptor}:{r['sequence']}\n" for r in candidates))
    config.update(state="candidates_and_msas_ready", n_sequences=len(candidates),
                  msa_records=msa_stats, msa_hashes=msa_hashes,
                  generation_seconds=time.perf_counter() - started,
                  software={"torch": torch.__version__, "numpy": np.__version__},
                  candidate_file_sha256=sha256(OUT / "candidates.csv"))
    dump_json(OUT / "protocol.json", config)
    print(f"Prepared {len(candidates)} sequences; elapsed {config['generation_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    generate()
