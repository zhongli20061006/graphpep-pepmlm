#!/usr/bin/env python3
"""Multi-parent MDM2 feedback panel: native vs predicted vs random site choice."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
import torch
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SVDSuperimposer import SVDSuperimposer
from transformers import AutoModelForMaskedLM, AutoTokenizer

import feedback_pilot as common
import step1_generate
import step3_graphpep as gp
from shared_receptor_msa import make_shared_msa

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/feedback_panel_mdm2_20261005"
PRIOR = ROOT / "output/experiment0/expanded_validation"
MODEL = Path("/home/g/pepmlm_pipeline/models/PepMLM-650M")
SEED = 20261009
TARGET = "ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVV"
PARENTS = [
    {"peptide_id": "P53WT_1YCR", "sequence": "ETFSDLWKLLPEN", "pdb": "1YCR",
     "full_sequence": "ETFSDLWKLLPEN", "native_path": PRIOR / "POSCTRL_MDM2_P53_native.pdb",
     "chain_rec": "A", "chain_pep": "B", "receptor_range": "25-109", "peptide_range": "17-29"},
    {"peptide_id": "PMI_3EQS", "sequence": "TSFAEYWNLLSP", "pdb": "3EQS",
     "full_sequence": "TSFAEYWNLLSP", "cif": OUT / "pdb/3EQS.cif", "chain_rec": "A", "chain_pep": "B",
     "receptor_range": "experimental construct; first receptor residue unresolved", "peptide_range": "12-mer; C-terminal P unresolved"},
    {"peptide_id": "P53LIKE_3G03", "sequence": "LTFEHYWAQLTS", "pdb": "3G03",
     "full_sequence": "LTFEHYWAQLTS", "cif": OUT / "pdb/3G03.cif", "chain_rec": "A", "chain_pep": "B",
     "receptor_range": "MDM2 canonical residues 25-109", "peptide_range": "12-mer; C-terminal S unresolved"},
]
AMINO_ACIDS = common.AMINO_ACIDS
PARSER = PDBParser(QUIET=True)
COMPARISON_PROTOCOL = {
    "model_gate": {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10., "minimum_distance_min": 1.5},
    "minimum_common_passing_models": 2, "minimum_score_gain": .05,
    "minimum_matched_passing_models_improved": 2,
    "afm_model_type": "alphafold2_multimer_v3", "afm_models": [1, 2, 3],
    "afm_recycles": 3, "afm_recycle_early_stop_tolerance": 0,
    "max_seq": 256, "max_extra_seq": 512, "templates": False,
    "relaxation": {"method": "ColabFold Amber restrained minimization", "max_iterations": 200,
                   "max_outer_iterations": 3, "tolerance": 2.39, "stiffness": 10., "use_gpu": False},
    "afm_seed": SEED,
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def standard_residues(chain) -> tuple[list, str]:
    residues, letters = [], []
    for residue in chain:
        if is_aa(residue, standard=True) and residue.resname in gp.THREE_TO_ONE:
            residues.append(residue)
            letters.append(gp.THREE_TO_ONE[residue.resname])
    return residues, "".join(letters)


def write_chain(residues: list, chain_id: str, starting_serial: int) -> tuple[str, int]:
    lines, seen, sequence, serial = [], set(), [], starting_serial
    index = 0
    last_resname = None
    for residue in residues:
        resname = residue.resname.upper()
        if resname not in gp.THREE_TO_ONE:
            continue
        atom_records = []
        for atom in residue:
            if atom.element.upper() in {"H", "D"} or atom.altloc not in {" ", "A"}:
                continue
            name = atom.name.strip()
            if name in seen:
                continue
            seen.add(name)
            atom_records.append(atom)
        if not atom_records:
            seen.clear()
            continue
        index += 1
        last_resname = resname
        sequence.append(gp.THREE_TO_ONE[resname])
        for atom in atom_records:
            x, y, z = atom.coord
            element = atom.element.strip().upper() or atom.name[0]
            atomname = atom.name.strip()
            atom_field = f" {atomname:<3s}" if len(atomname) < 4 and len(element) == 1 else f"{atomname:>4s}"
            lines.append(f"ATOM  {serial:5d} {atom_field} {resname:>3s} {chain_id}{index:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{atom.occupancy or 1.:6.2f}{atom.bfactor:6.2f}"
                         f"          {element:>2s}\n")
            serial += 1
        seen.clear()
    if not index or last_resname is None:
        raise ValueError(f"no standard residues to write for chain {chain_id}")
    lines.append(f"TER   {serial:5d}      {last_resname:>3s} {chain_id}{index:4d}\n")
    return "".join(lines), index


def model_atoms(lines: list[str], chain: str) -> tuple[str, list, np.ndarray, np.ndarray]:
    residues, seen = [], set()
    for line in lines:
        if not line.startswith("ATOM") or line[21] != chain:
            continue
        key = (line[22:26], line[26])
        if key not in seen:
            seen.add(key)
            residues.append(line[17:20].strip())
    seq = "".join(gp.THREE_TO_ONE[r] for r in residues)
    xyz, positions = [], []
    for line in lines:
        if not line.startswith("ATOM") or line[21] != chain or line[76:78].strip().upper() in {"H", "D"}:
            continue
        xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        positions.append(int(line[22:26]))
    return seq, residues, np.asarray(xyz, dtype=float), np.asarray(positions, dtype=int)


def extract_native(parent: dict, parser: MMCIFParser) -> dict:
    path = OUT / "native_complexes" / f"{parent['peptide_id']}_native.pdb"
    path.parent.mkdir(parents=True, exist_ok=True)
    if parent["pdb"] == "1YCR":
        source = parent["native_path"]
        if digest(source) != json.loads((PRIOR / "positive_control.json").read_text())["native_complex_sha256"]:
            raise ValueError("1YCR native parent changed")
        native_structure = PARSER.get_structure("1YCR", str(source))
        rec_res, rec_seq = standard_residues(native_structure[0][parent["chain_rec"]])
        pep_res, pep_seq = standard_residues(native_structure[0][parent["chain_pep"]])
        if pep_seq != parent["sequence"] or rec_seq != TARGET:
            raise ValueError("1YCR parent chains differ from the frozen sequence")
        path.write_bytes(source.read_bytes())
        result = {key: value for key, value in parent.items() if key not in {"native_path", "cif"}}
        result.update(source_coordinate_path=str(source), native_path=str(path),
                      native_receptor_sequence=rec_seq, native_peptide_coordinate_sequence=pep_seq,
                      native_peptide_unresolved_positions=[], source_sha256=digest(source))
        return result
    structure = parser.get_structure(parent["pdb"], str(parent["cif"]))
    chains = structure[0]
    receptor_residues, rec_seq = standard_residues(chains[parent["chain_rec"]])
    peptide_residues, pep_seq = standard_residues(chains[parent["chain_pep"]])
    if parent["pdb"] == "3EQS":
        if TARGET[1:] != rec_seq:
            raise ValueError(f"3EQS resolved MDM2 chain unexpected: {rec_seq}")
    else:
        if not rec_seq.startswith(TARGET) or len(rec_seq) < len(TARGET):
            raise ValueError(f"3G03 resolved MDM2 chain unexpected: {rec_seq}")
        receptor_residues = receptor_residues[:len(TARGET)]
        rec_seq = TARGET
    sequence = parent["sequence"]
    if not sequence.startswith(pep_seq) or len(pep_seq) < 10:
        raise ValueError(f"PDB peptide coordinates {pep_seq} do not match documented {sequence}")
    rec_lines, rec_count = write_chain(receptor_residues, "A", 1)
    pep_lines, pep_count = write_chain(peptide_residues, "B", 1 + sum(line.startswith("ATOM") for line in rec_lines))
    if rec_count != len(rec_seq) or pep_count != len(pep_seq):
        raise ValueError("Unexpected chain residue loss while writing normalized PDB")
    path.write_text(rec_lines + pep_lines + "END\n")
    result = {key: value for key, value in parent.items() if key not in {"native_path", "cif"}}
    result.update(source_coordinate_path=str(parent["cif"]), native_path=str(path),
                  native_receptor_sequence=rec_seq, native_peptide_coordinate_sequence=pep_seq,
                  native_peptide_unresolved_positions=list(range(len(pep_seq) + 1, len(sequence) + 1)),
                  source_sha256=digest(parent["cif"]))
    return result


def run_graphpep(manifest: pd.DataFrame, directory: Path) -> pd.DataFrame:
    pairs = gp.build_inputs(manifest, str(directory))
    results = []
    mapper = gp.ROOT / "pipeline/graphpep_patches/resmap_compat.py"
    signature_paths = [gp.GP / name for name in ["GraphPep.sh", "bin/core.py", "bin/data.py", "bin/predict.py"]]
    signature_paths += [mapper, Path(gp.__file__)]
    signature = hashlib.sha256("".join(digest(path) for path in signature_paths).encode()).hexdigest()
    score_path = directory / "scores.csv"
    completed = pd.read_csv(score_path).to_dict("records") if score_path.is_file() else []
    cache = {(str(row["peptide_id"]), str(row["model_id"]), str(row["complex_path"])): row for row in completed}
    for index, pair in enumerate(pairs, 1):
        row = pair["row"]
        key = (str(row["peptide_id"]), str(row["model_id"]), str(row["complex_path"]))
        structure_hash = digest(Path(row["complex_path"]))
        saved = cache.get(key, {})
        reusable = (saved.get("graphpep_status") == "ok" and
                    saved.get("structure_sha256") == structure_hash and
                    saved.get("graphpep_source_signature") == signature and
                    Path(str(saved.get("contact_feedback_path", ""))).is_file())
        print(f"GraphPep {index}/{len(pairs)}{' cached' if reusable else ''}: {row['peptide_id']} / {row['model_id']}", flush=True)
        result = {**saved, **row} if reusable else gp.score_pair(pair, gp.GP, mapper, .7)
        result.update(structure_sha256=structure_hash, graphpep_source_signature=signature)
        results.append(result)
        pd.DataFrame(results).to_csv(score_path, index=False)
    return pd.DataFrame(results)


def make_parent_inputs() -> None:
    if (OUT / "parent_metadata.json").exists():
        raise FileExistsError("Frozen parent inputs already exist")
    parser = MMCIFParser(QUIET=True)
    parents = [extract_native(parent, parser) for parent in PARENTS]
    cached = (PRIOR / "results/POSCTRL_MDM2_P53.a3m").read_text()
    msa_folder = OUT / "parent_msa"
    msa_folder.mkdir(parents=True, exist_ok=True)
    fasta_lines, native_manifest = [], []
    for parent in parents:
        msa, stats = make_shared_msa(cached, TARGET, PARENTS[0]["sequence"], parent["sequence"])
        msa_path = msa_folder / f"{parent['peptide_id']}.a3m"
        msa_path.write_text(msa)
        parent["msa_path"] = str(msa_path)
        parent["msa_sha256"] = digest(msa_path)
        parent["msa_stats"] = stats
        fasta_lines += [f">{parent['peptide_id']}\n", f"{TARGET}:{parent['sequence']}\n"]
        native_complex = Path(parent["native_path"])
        native_manifest.append({"peptide_id": parent["peptide_id"],
                                "peptide_sequence": parent["native_peptide_coordinate_sequence"],
                                "receptor_sequence": parent["native_receptor_sequence"],
                                "source": "experimental_native", "model_id": parent["pdb"],
                                "complex_path": str(native_complex)})
    (OUT / "parent_complexes.fasta").write_text("".join(fasta_lines))
    if len({parent["msa_stats"]["receptor_alignment_sha256"] for parent in parents}) != 1:
        raise ValueError("Receptor alignments differ between the three parent MSAs")
    pd.DataFrame(native_manifest).to_csv(OUT / "native_manifest.csv", index=False)
    native_scores = run_graphpep(pd.DataFrame(native_manifest), OUT / "native_graphpep")
    native_scores.to_csv(OUT / "native_scores.csv", index=False)
    for parent in parents:
        row = native_scores[native_scores.peptide_id == parent["peptide_id"]].iloc[0]
        parent["native_feedback_path"] = row.residue_feedback_path
        parent["native_contact_feedback_path"] = row.contact_feedback_path
        parent["native_graphpep_score"] = float(row.binding_score)
        parent["native_contacts"] = int(row.n_contacts)
    write_json(OUT / "parent_metadata.json", {"stage": "native_controls_ready", "seed": SEED,
                                               "receptor_sequence": TARGET, "parents": parents,
                                               "native_scores_sha256": digest(OUT / "native_scores.csv"),
                                               "comparison_protocol": {
                                                   "model_gate": {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10., "minimum_distance_min": 1.5},
                                                   "minimum_common_passing_models": 2, "minimum_score_gain": .05,
                                                   "minimum_matched_passing_models_improved": 2,
                                                   "afm_model_type": "alphafold2_multimer_v3", "afm_models": [1, 2, 3],
                                                   "afm_recycles": 3, "afm_recycle_early_stop_tolerance": 0,
                                                   "max_seq": 256, "max_extra_seq": 512, "templates": False,
                                                   "relaxation": False, "afm_seed": SEED},
                                               "source_a3m_sha256": digest(PRIOR / "results/POSCTRL_MDM2_P53.a3m"),
                                               "native_graphpep_source_hashes": {
                                                   str(gp.GP / name): digest(gp.GP / name)
                                                   for name in ["GraphPep.sh", "bin/core.py", "bin/data.py"]}})
    print(f"Prepared {len(parents)} native parents and shared-receptor MSAs", flush=True)


def parent_rows() -> list[dict]:
    return json.loads((OUT / "parent_metadata.json").read_text())["parents"]


def predictions(folder: Path, parent: dict) -> list[Path]:
    state = "relaxed" if "relaxed" in folder.name else "unrelaxed"
    found = sorted(folder.glob(f"{parent['peptide_id']}_{state}_rank_*_model_*_seed_{SEED}.pdb"))
    if len(found) != 3:
        raise ValueError(f"Expected 3 AFM structures for {parent['peptide_id']}, found {len(found)}")
    return found


def confidence_file(path: Path) -> Path:
    return path.with_name(re.sub(r"_(unrelaxed|relaxed)_", "_scores_", path.name)).with_suffix(".json")


def profile(path: str, length: int) -> list[dict]:
    frame = pd.read_csv(path)
    if "peptide_position" not in frame:
        frame = frame.rename(columns={frame.columns[0]: "peptide_position"})
    required = {"peptide_position", "count", "mean"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Unexpected GraphPep feedback format: {frame.columns.tolist()}")
    return [{"peptide_position": i, "count": int(frame.loc[frame.peptide_position == i, "count"].iloc[0])
             if (frame.peptide_position == i).any() else 0,
             "mean": float(frame.loc[frame.peptide_position == i, "mean"].iloc[0])
             if (frame.peptide_position == i).any() and pd.notna(frame.loc[frame.peptide_position == i, "mean"].iloc[0]) else ""}
            for i in range(1, length + 1)]


def parent_structure_scores() -> None:
    metadata = json.loads((OUT / "parent_metadata.json").read_text())
    manifest = []
    for parent in metadata["parents"]:
        for index, path in enumerate(predictions(OUT / "parent_relaxed", parent), 1):
            number = int(re.search(r"_model_(\d+)_seed_", path.name).group(1))
            manifest.append({"peptide_id": parent["peptide_id"], "peptide_sequence": parent["sequence"],
                             "receptor_sequence": TARGET, "source": "predicted_parent", "model_id": f"model_{number}",
                             "complex_path": str(path)})
    scores = run_graphpep(pd.DataFrame(manifest), OUT / "parent_relaxed_graphpep")
    score_path = OUT / "parent_model_scores.csv"
    scores.to_csv(score_path, index=False)
    model_data = []
    for row in scores.to_dict("records"):
        pred = Path(row["complex_path"])
        confidence = confidence_file(pred)
        values = json.loads(confidence.read_text())
        nrec = len(TARGET)
        rec = residue_atom_data(pred, "A")
        pep = residue_atom_data(pred, "B")
        pairs6, minimum, clashes = distances(rec, pep, 6.)
        contacting_rec = sorted({r - 1 for r, _ in pairs6})
        pae = np.asarray(values["pae"])
        peptide_indices = list(range(nrec, nrec + len(row["peptide_sequence"])))
        interface_pae = float(np.median(np.concatenate([
            pae[np.ix_(contacting_rec, peptide_indices)].ravel(),
            pae[np.ix_(peptide_indices, contacting_rec)].ravel()])) if contacting_rec else float("nan"))
        peptide_plddt = float(np.mean(values["plddt"][nrec:]))
        limits = {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10., "minimum_distance_min": 1.5}
        reasons = []
        if float(values["iptm"]) < limits["iptm_min"]: reasons.append("low_iptm")
        if peptide_plddt < limits["peptide_plddt_min"]: reasons.append("low_peptide_plddt")
        if not math.isfinite(interface_pae) or interface_pae > limits["interface_pae_max"]: reasons.append("high_or_missing_interface_pae")
        if minimum < limits["minimum_distance_min"]: reasons.append("heavyatom_clash")
        if not pairs6: reasons.append("no_interface_contacts")
        model_data.append({"peptide_id": row["peptide_id"], "model_id": row["model_id"],
                           "sequence": row["peptide_sequence"], "iptm": values["iptm"],
                           "peptide_plddt": peptide_plddt, "interface_pae": interface_pae,
                           "contact_pairs_6A": len(pairs6), "minimum_heavyatom_distance": minimum,
                           "clashing_atom_pairs_lt_1p5A": clashes,
                           "binding_score": row["binding_score"], "mean_confidence": row["mean_conf"],
                           "contact_count": row["n_contacts"], "structure_path": str(pred),
                           "residue_feedback_path": row["residue_feedback_path"],
                           "structure_pass": not reasons, "quality_fail_reasons": ";".join(reasons)})
    pd.DataFrame(model_data).to_csv(OUT / "parent_model_metrics.csv", index=False)
    for parent in metadata["parents"]:
        group = [row for row in model_data if row["peptide_id"] == parent["peptide_id"]]
        contact_moments = np.zeros(len(parent["sequence"]), dtype=float)
        total_contacts = np.zeros(len(parent["sequence"]), dtype=int)
        usable_group = [row for row in group if row["structure_pass"]]
        if len(usable_group) < 2:
            raise ValueError(f"Only {len(usable_group)}/3 parent AFM structures pass quality gates for {parent['peptide_id']}; predicted-feedback arm is not interpretable")
        for row in usable_group:
            for point in profile(row["residue_feedback_path"], len(parent["sequence"])):
                pos = point["peptide_position"] - 1
                if point["count"]:
                    contact_moments[pos] += point["count"] * float(point["mean"])
                    total_contacts[pos] += point["count"]
        predicted_profile = [{"peptide_position": pos + 1, "count": int(count),
                              "mean": float(contact_moments[pos] / count) if count else ""}
                             for pos, count in enumerate(total_contacts)]
        pred_path = OUT / f"{parent['peptide_id']}_predicted_ensemble_feedback.csv"
        pd.DataFrame(predicted_profile).to_csv(pred_path, index=False)
        native_profile = profile(parent["native_feedback_path"], len(parent["sequence"]))
        native_eligible, native_sites = common.select_positions(native_profile, 2)
        pred_eligible, pred_sites = common.select_positions(predicted_profile, 2)
        comparison = pd.DataFrame(native_profile).rename(columns={"count": "native_contacts", "mean": "native_mean"})
        prediction_frame = pd.DataFrame(predicted_profile).rename(columns={"count": "predicted_contacts", "mean": "predicted_mean"})
        comparison = comparison.merge(prediction_frame, on="peptide_position", how="outer")
        comparison.to_csv(OUT / f"{parent['peptide_id']}_feedback_comparison.csv", index=False)
        paired = comparison[(comparison.native_contacts > 0) & (comparison.predicted_contacts > 0)]
        native_values = pd.to_numeric(paired.native_mean)
        predicted_values = pd.to_numeric(paired.predicted_mean)
        pearson = native_values.corr(predicted_values)
        spearman = native_values.rank().corr(predicted_values.rank())
        parent.update(native_profile_path=parent["native_feedback_path"],
                      predicted_ensemble_profile_path=str(pred_path),
                      native_eligible_positions=native_eligible, native_guided_positions=sorted(native_sites),
                      predicted_eligible_positions=pred_eligible, predicted_guided_positions=sorted(pred_sites),
                      predicted_profile_coverage=int(sum(p["count"] > 0 for p in predicted_profile)),
                      predicted_feedback_model_ids=[row["model_id"] for row in usable_group],
                      native_predicted_feedback_common_positions=len(paired),
                      native_predicted_feedback_pearson=float(pearson) if pd.notna(pearson) else None,
                      native_predicted_feedback_spearman=float(spearman) if pd.notna(spearman) else None,
                      selected_site_overlap=len(set(native_sites) & set(pred_sites)))
    metadata.update(stage="parent_structures_scored_feedback_frozen", parent_model_score_sha256=digest(score_path),
                    parents=metadata["parents"],
                    sampling_protocol={"groups": ["native_guided", "native_random", "predicted_guided", "predicted_random"],
                                       "children_per_group_parent": 3, "mutations_per_child": 2,
                                       "top_k": 3, "temperature": 1.0,
                                       "random_sites_sampled_from_each_feedback_source_eligible_set": True,
                                       "generation_seed": SEED,
                                       "predicted_feedback": "contact-weighted mean across quality-passing AFM parent models 1-3"})
    write_json(OUT / "parent_metadata.json", metadata)


def draw_child(parent: str, positions: list[int], model, tokenizer, device: str,
               seed: int, seen: set[str], top_k: int = 3) -> tuple[str, dict]:
    if top_k < 2:
        raise ValueError("top_k must leave multiple amino-acid choices at every masked site")
    masked = "".join("<mask>" if index in positions else residue for index, residue in enumerate(parent, 1))
    encoded = tokenizer(TARGET + masked, return_tensors="pt").to(device)
    mask_positions = (encoded.input_ids[0] == tokenizer.mask_token_id).nonzero(as_tuple=True)[0].tolist()
    expected = [len(TARGET) + pos for pos in positions]
    if mask_positions != expected:
        raise ValueError(f"PepMLM mask offset mismatch: {mask_positions} != {expected}")
    canonical = [tokenizer.convert_tokens_to_ids(aa) for aa in AMINO_ACIDS]
    with torch.inference_mode():
        logits = model(**encoded).logits[0, mask_positions].float().cpu()
    choices = []
    for index, pos in enumerate(positions):
        scores = logits[index, canonical].clone()
        scores[AMINO_ACIDS.index(parent[pos - 1])] = -torch.inf
        top, ids = scores.topk(top_k)
        probs = torch.softmax(top, 0).tolist()
        choices.append([(AMINO_ACIDS[aa_id], prob) for aa_id, prob in zip(ids.tolist(), probs)])
    alternatives = []
    for choice in itertools.product(*choices):
        letters = list(parent)
        for pos, (aa, _) in zip(positions, choice):
            letters[pos - 1] = aa
        sequence = "".join(letters)
        if sequence not in seen:
            alternatives.append((sequence, math.prod(prob for _, prob in choice)))
    probs = np.asarray([prob for _, prob in alternatives], dtype=float)
    probs /= probs.sum()
    rng = np.random.default_rng(seed)
    selected = int(rng.choice(len(alternatives), p=probs))
    child = alternatives[selected][0]
    common.validate_child(parent, child, positions)
    seen.add(child)
    return child, {"seed": seed, "positions": positions, "alternatives": alternatives,
                   "normalized_probabilities": probs.tolist(), "selected_index": selected}


def generate_children() -> None:
    metadata = json.loads((OUT / "parent_metadata.json").read_text())
    if metadata["stage"] != "parent_structures_scored_feedback_frozen":
        raise ValueError("Parent structure feedback and site choices must be frozen first")
    outcsv = OUT / "candidates.csv"
    if outcsv.exists():
        raise FileExistsError("Candidate set is already frozen")
    metadata["comparison_protocol"] = COMPARISON_PROTOCOL
    verify_parent_metrics(metadata)
    snapshot = OUT / "source_snapshot"
    snapshot.mkdir(exist_ok=True)
    for source in [Path(__file__), ROOT / "pipeline/shared_receptor_msa.py", ROOT / "pipeline/step1_generate.py", ROOT / "pipeline/relax_feedback_panel.py"]:
        (snapshot / source.name).write_bytes(source.read_bytes())
    write_json(OUT / "protocol.json", {
        "seed": SEED, "comparison": COMPARISON_PROTOCOL,
        "sampling": metadata["sampling_protocol"],
        "parents": [{key: parent[key] for key in ["peptide_id", "sequence", "pdb", "native_guided_positions",
                     "predicted_guided_positions", "native_eligible_positions", "predicted_eligible_positions",
                     "predicted_feedback_model_ids"]} for parent in metadata["parents"]],
        "parent_model_scores_sha256": metadata["parent_model_score_sha256"],
        "source_hashes": {str(path): digest(path) for path in [Path(__file__), ROOT / "pipeline/shared_receptor_msa.py", ROOT / "pipeline/step1_generate.py"]},
        "source_snapshot_hashes": {str(path): digest(path) for path in snapshot.glob("*.py")},
    })
    torch.manual_seed(SEED)
    tokenizer = AutoTokenizer.from_pretrained(MODEL, local_files_only=True)
    model = AutoModelForMaskedLM.from_pretrained(MODEL, torch_dtype=torch.float16, local_files_only=True).eval()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(device)
    rows, sampling_log, fasta = [], [], []
    for parent_index, parent_row in enumerate(metadata["parents"]):
        parent = parent_row["sequence"]
        fasta += [f">{parent_row['peptide_id']}\n", f"{TARGET}:{parent}\n"]
        rows.append({"peptide_id": parent_row["peptide_id"], "parent_id": parent_row["peptide_id"],
                     "arm": "parent", "replicate": 0, "sequence": parent, "length": len(parent),
                     "mask_positions": "", "mutations": "", "mutation_count": 0, "ppl_score": np.nan})
        rng = np.random.default_rng(SEED + parent_index * 100)
        groups = [
            ("NATIVE_GUIDED", parent_row["native_guided_positions"]),
            ("NATIVE_RANDOM", None),
            ("PREDICTED_GUIDED", parent_row["predicted_guided_positions"]),
            ("PREDICTED_RANDOM", None),
        ]
        for arm_index, (arm, guided) in enumerate(groups):
            eligible = (parent_row["native_eligible_positions"] if arm.startswith("NATIVE")
                        else parent_row["predicted_eligible_positions"])
            seen: set[str] = set()
            for replicate in range(1, 4):
                positions = sorted(guided) if guided is not None else sorted(rng.choice(eligible, 2, replace=False).tolist())
                sample_seed = SEED + parent_index * 1000 + arm_index * 10 + replicate
                child, trace = draw_child(parent, positions, model, tokenizer, device, sample_seed, seen)
                mutations = ";".join(f"{parent[pos-1]}{pos}{child[pos-1]}" for pos in positions)
                ident = f"{arm}_{parent_row['peptide_id']}_{replicate:02d}"
                rows.append({"peptide_id": ident, "parent_id": parent_row["peptide_id"], "arm": arm.lower(),
                             "replicate": replicate, "sequence": child, "length": len(child),
                             "mask_positions": ";".join(map(str, positions)), "mutations": mutations,
                             "mutation_count": len(positions), "ppl_score": np.nan, "sampling_seed": sample_seed})
                sampling_log.append({"peptide_id": ident, "parent_id": parent_row["peptide_id"], "arm": arm.lower(),
                                     "guided_positions": guided, "eligible_positions": eligible,
                                     "selected_positions": positions, "mutation_summary": mutations,
                                     "sampling": trace, "sequence": child})
                fasta += [f">{ident}\n", f"{TARGET}:{child}\n"]
                print(f"{ident}: {child} ({mutations})", flush=True)
    frame = pd.DataFrame(rows)
    device_name = torch.cuda.get_device_name(0) if device == "cuda" else "CPU"
    for parent in metadata["parents"]:
        indices = frame.index[frame.parent_id == parent["peptide_id"]].tolist()
        sequences = frame.loc[indices, "sequence"].tolist()
        target_ids = torch.tensor(tokenizer(TARGET)["input_ids"][:-1])
        values = step1_generate.pseudo_ppl_chunk(model, tokenizer, target_ids, sequences, device, len(parent["sequence"]), 16)
        frame.loc[indices, "ppl_score"] = values
    frame.to_csv(outcsv, index=False)
    (OUT / "candidate_complexes.fasta").write_text("".join(fasta))
    write_json(OUT / "sampling_log.json", {"seed": SEED, "models": [1, 2, 3], "proposals": sampling_log})
    msa_dir = OUT / "candidate_msa"
    msa_dir.mkdir(exist_ok=True)
    offspring_msa_dir = OUT / "offspring_msa"
    offspring_msa_dir.mkdir(exist_ok=True)
    cached = (PRIOR / "results/POSCTRL_MDM2_P53.a3m").read_text()
    msa_hashes, msa_records = {}, []
    for row in rows:
        msa, stats = make_shared_msa(cached, TARGET, PARENTS[0]["sequence"], row["sequence"])
        path = msa_dir / f"{row['peptide_id']}.a3m"
        path.write_text(msa)
        msa_hashes[str(path)] = digest(path)
        msa_records.append({"peptide_id": row["peptide_id"], **stats})
        if row["arm"] != "parent":
            (offspring_msa_dir / path.name).write_text(msa)
    if len({row["receptor_alignment_sha256"] for row in msa_records}) != 1:
        raise ValueError("Receptor alignments differ across candidate inputs")
    if len(set(msa_hashes.values())) == 1:
        raise ValueError("Every MSA is unexpectedly byte-identical despite distinct peptide queries")
    metadata.update(stage="offspring_and_msas_ready", candidate_sha256=digest(outcsv),
                    protocol_sha256=digest(OUT / "protocol.json"),
                    candidate_msa_hashes=msa_hashes, msa_records=msa_records, pepmlm_device=device_name,
                    source_hashes={str(Path(row)): digest(Path(row)) for row in metadata["native_graphpep_source_hashes"]},
                    software={"torch": torch.__version__, "numpy": np.__version__})
    write_json(OUT / "parent_metadata.json", metadata)


def distances(rec, pep, cutoff: float = 6.0) -> tuple[set[tuple[int, int]], float, int]:
    rec_xyz, rec_pos = rec[2], rec[3]
    pep_xyz, pep_pos = pep[2], pep[3]
    matrix = np.sqrt(np.sum((rec_xyz[:, None] - pep_xyz[None, :]) ** 2, axis=2))
    a, b = np.where(matrix <= cutoff)
    pairs = {(int(rec_pos[x]), int(pep_pos[y])) for x, y in zip(a, b)}
    return pairs, float(matrix.min()), int((matrix < 1.5).sum())


def residue_atom_data(path: Path, chain: str) -> tuple[str, list, np.ndarray, np.ndarray, dict[int, np.ndarray]]:
    structure = PARSER.get_structure(path.stem, str(path))
    residues, sequence = standard_residues(structure[0][chain])
    xyz, positions, cas = [], [], {}
    for pos, residue in enumerate(residues, 1):
        if residue.has_id("CA"):
            cas[pos] = residue["CA"].coord.copy()
        for atom in residue:
            if atom.element.upper() not in {"H", "D"}:
                xyz.append(atom.coord)
                positions.append(pos)
    return sequence, residues, np.asarray(xyz), np.asarray(positions), cas


def sequence_map(reference: str, observed: str) -> dict[int, int]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2.0, -3.0
    aligner.open_gap_score, aligner.extend_gap_score = -4.0, -0.5
    mapping = {}
    alignment = aligner.align(reference, observed)[0]
    left, right = alignment.aligned
    for (a0, a1), (b0, b1) in zip(left, right):
        for a, b in zip(range(a0, a1), range(b0, b1)):
            if reference[a] == observed[b]:
                mapping[b + 1] = a + 1
    return mapping


def metrics(parent: dict, row: dict, path: Path, native_path: Path, native_contact_positions: set[tuple[int, int]]) -> dict:
    rec = residue_atom_data(path, "A")
    pep = residue_atom_data(path, "B")
    if rec[0] != TARGET or pep[0] != row["sequence"]:
        raise ValueError(f"Predicted chain sequence mismatch in {path}")
    nr, np_ = len(TARGET), len(row["sequence"])
    conf_path = confidence_file(path)
    conf = json.loads(conf_path.read_text())
    pae, plddt = np.asarray(conf["pae"]), np.asarray(conf["plddt"])
    if pae.shape != (nr + np_, nr + np_) or plddt.shape != (nr + np_,):
        raise ValueError("Predicted confidence shape disagrees with chain lengths")
    pairs5, minimum, clashes = distances(rec, pep, 5.0)
    pairs6, _, _ = distances(rec, pep, 6.0)
    contacting_rec = sorted({r - 1 for r, _ in pairs6})
    pep_index = list(range(nr, nr + np_))
    interface_pae = float(np.median(np.concatenate([pae[np.ix_(contacting_rec, pep_index)].ravel(),
                                                     pae[np.ix_(pep_index, contacting_rec)].ravel()]))) if contacting_rec else float("nan")
    native = residue_atom_data(native_path, "A")
    native_pep = residue_atom_data(native_path, "B")
    rec_map = sequence_map(rec[0], native[0])  # native position -> canonical AFM receptor position
    sup = SVDSuperimposer()
    fixed, moving = [], []
    for native_pos, model_pos in rec_map.items():
        if native_pos in native[4] and model_pos in rec[4]:
            fixed.append(native[4][native_pos])
            moving.append(rec[4][model_pos])
    if len(fixed) < 20:
        raise ValueError(f"Too few aligned MDM2 residues in {native_path.name}: {len(fixed)}")
    sup.set(np.asarray(fixed), np.asarray(moving))
    sup.run()
    rotation, translation = sup.get_rotran()
    common_positions = min(len(native_pep[4]), len(pep[4]))
    aligned = np.asarray([pep[4][i] @ rotation + translation for i in range(1, common_positions + 1)])
    crystal = np.asarray([native_pep[4][i] for i in range(1, common_positions + 1)])
    peptide_rmsd = float(np.sqrt(np.mean(np.sum((aligned - crystal) ** 2, axis=1)))) if common_positions else float("nan")
    mapped_native = {(rec_map[r], p) for r, p in native_contact_positions if r in rec_map}
    native_recall = len(pairs5 & mapped_native) / len(mapped_native) if mapped_native else float("nan")
    limits, reasons = parent["gate"], []
    iptm = float(conf["iptm"])
    pep_plddt = float(plddt[nr:].mean())
    if iptm < limits["iptm_min"]: reasons.append("low_iptm")
    if pep_plddt < limits["peptide_plddt_min"]: reasons.append("low_peptide_plddt")
    if not math.isfinite(interface_pae) or interface_pae > limits["interface_pae_max"]: reasons.append("high_or_missing_interface_pae")
    if minimum < limits["minimum_distance_min"]: reasons.append("heavyatom_clash")
    if not pairs6: reasons.append("no_interface_contacts")
    return {**row, "peptide_sequence": row["sequence"], "receptor_sequence": TARGET,
            "source": "predicted_parent" if row["arm"] == "parent" else "predicted_offspring",
            "model_id": path.name, "model_path": str(path), "confidence_path": str(conf_path),
            "iptm": iptm, "ptm": float(conf["ptm"]), "peptide_plddt": pep_plddt,
            "interface_pae": interface_pae, "contact_pairs_5A": len(pairs5), "contact_pairs_6A": len(pairs6),
            "minimum_heavyatom_distance": minimum, "clashing_atom_pairs_lt_1p5A": clashes,
            "contact_model_jaccard": np.nan, "native_contact_recall": native_recall,
            "receptor_CA_RMSD_after_alignment_A": float(sup.get_rms()),
            "peptide_CA_RMSD_after_receptor_alignment_A": peptide_rmsd,
            "structure_pass": not reasons, "quality_fail_reasons": ";".join(reasons),
            "model_number": int(re.search(r"_model_(\d+)_seed_", path.name).group(1))}


def score_bounds(row: dict) -> None:
    if row["graphpep_status"] == "no_contact":
        if row["contact_pairs_6A"]:
            raise ValueError("GraphPep no-contact status disagrees with geometry")
        return
    values = pd.read_csv(row["contact_feedback_path"])["CONFIDENT_SCORES"].to_numpy()
    low = sum(w * int((values - .00005 > cut).sum()) for w, cut in [(1, .5), (2, .7), (4, .9)]) / 7
    high = sum(w * int((values + .00005 > cut).sum()) for w, cut in [(1, .5), (2, .7), (4, .9)]) / 7
    expected = (-math.log1p(high), -math.log1p(low))
    if not expected[0] - .000051 <= float(row["binding_score"]) <= expected[1] + .000051:
        raise ValueError("GraphPep raw score disagrees with independent formula check")
    if len(values) != int(row["n_contacts"]) or len(values) != int(row["contact_pairs_6A"]):
        raise ValueError("GraphPep output contact count disagrees with geometry")


def verify_parent_metrics(metadata: dict) -> None:
    """Cross-check nine parent geometries and their frozen evaluator outputs."""
    saved_scores = pd.read_csv(OUT / "parent_model_scores.csv")
    rows = []
    for parent in metadata["parents"]:
        native = Path(parent["native_path"])
        native_contacts = collect_native_contact_positions(parent)
        receptor_mapping = sequence_map(TARGET, parent["native_receptor_sequence"])
        offset = 1 if parent["pdb"] == "3EQS" else 0
        expected_mapping = {i: i + offset for i in range(1, len(parent["native_receptor_sequence"]) + 1)}
        if receptor_mapping != expected_mapping:
            raise ValueError("Experimental receptor positions do not map to the expected MDM2 construct")
        spec = {**parent, "gate": COMPARISON_PROTOCOL["model_gate"]}
        for path in predictions(OUT / "parent_relaxed", parent):
            row = metrics(spec, {"peptide_id": parent["peptide_id"], "parent_id": parent["peptide_id"],
                                 "sequence": parent["sequence"], "arm": "parent"}, path, native, native_contacts)
            saved = saved_scores[saved_scores.complex_path == str(path)]
            if len(saved) != 1:
                raise ValueError("Predicted parent score is missing or duplicated")
            row.update(saved.iloc[0].to_dict())
            score_bounds(row)
            rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "parent_structure_validation.csv", index=False)


def parent_specs_from_outputs(metadata: dict, candidates: pd.DataFrame) -> dict[str, dict]:
    specs = {}
    for parent in metadata["parents"]:
        gate = {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10., "minimum_distance_min": 1.5}
        specs[parent["peptide_id"]] = {**parent, "gate": gate}
    return specs


def analyze() -> None:
    (OUT / "analysis").mkdir(parents=True, exist_ok=True)
    metadata = json.loads((OUT / "parent_metadata.json").read_text())
    if metadata.get("stage") != "offspring_and_msas_ready":
        raise ValueError("Expected generated children and frozen candidate MSAs")
    if digest(OUT / "candidates.csv") != metadata["candidate_sha256"]:
        raise ValueError("Frozen candidate manifest changed")
    if digest(OUT / "parent_model_scores.csv") != metadata["parent_model_score_sha256"]:
        raise ValueError("Frozen predicted-parent scores changed")
    if digest(OUT / "protocol.json") != metadata["protocol_sha256"]:
        raise ValueError("Frozen comparison protocol changed")
    for parent in metadata["parents"]:
        if digest(Path(parent["source_coordinate_path"])) != parent["source_sha256"]:
            raise ValueError(f"Experimental source coordinates changed for {parent['peptide_id']}")
    for path, expected in metadata["candidate_msa_hashes"].items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Frozen MSA changed: {path}")
    candidates = pd.read_csv(OUT / "candidates.csv", keep_default_na=False)
    for parent in metadata["parents"]:
        for arm, group in candidates[(candidates.parent_id == parent["peptide_id"]) & (candidates.arm != "parent")].groupby("arm"):
            if len(group) != 3 or group.sequence.nunique() != 3:
                raise ValueError(f"Unexpected children or within-arm duplicates: {parent['peptide_id']} {arm}")
            for row in group.to_dict("records"):
                positions = [int(pos) for pos in str(row["mask_positions"]).split(";")]
                common.validate_child(parent["sequence"], row["sequence"], positions)
                if len(positions) != 2 or int(row["mutation_count"]) != 2:
                    raise ValueError("Mutation budget changed")
    specs = parent_specs_from_outputs(metadata, candidates)
    model_rows, contact_sets = [], {}
    for candidate in candidates.to_dict("records"):
        parent = specs[candidate["parent_id"]]
        ident = candidate["peptide_id"]
        folder = OUT / "parent_relaxed" if candidate["arm"] == "parent" else OUT / "relaxed_structures"
        found = predictions(folder, {**candidate, "peptide_id": ident})
        for path in found:
            result = metrics(parent, candidate, path, Path(parent["native_path"]),
                             {(int(a), int(b)) for a, b in parent["native_contact_positions"]})
            model_rows.append(result)
            pairs = distances(residue_atom_data(path, "A"), residue_atom_data(path, "B"), 5.)[0]
            contact_sets[(ident, result["model_number"])] = pairs
    for result in model_rows:
        pairs = contact_sets[(result["peptide_id"], result["model_number"])]
        similarities = []
        for number in (1, 2, 3):
            if number == result["model_number"]: continue
            other = contact_sets[(result["peptide_id"], number)]
            union = pairs | other
            similarities.append(len(pairs & other) / len(union) if union else 1.)
        result["contact_model_jaccard"] = float(np.mean(similarities))
    metrics_frame = pd.DataFrame(model_rows)
    metrics_frame.to_csv(OUT / "analysis/model_metrics.csv", index=False)
    manifest = metrics_frame.rename(columns={"model_path": "complex_path"}).copy()
    offspring_scores = run_graphpep(manifest[manifest.arm != "parent"], OUT / "analysis/offspring_graphpep")
    cached_parents = pd.read_csv(OUT / "parent_model_scores.csv")
    parent_scores = []
    for row in manifest[manifest.arm == "parent"].to_dict("records"):
        saved = cached_parents[cached_parents.complex_path == row["complex_path"]]
        if len(saved) != 1:
            raise ValueError("Expected one frozen parent score per structure")
        parent_scores.append({**saved.iloc[0].to_dict(), **row})
    score_frame = pd.concat([pd.DataFrame(parent_scores), offspring_scores], ignore_index=True)
    for row in score_frame.to_dict("records"):
        score_bounds(row)
    score_frame.to_csv(OUT / "analysis/model_results.csv", index=False)
    summaries = []
    for ident, group in score_frame[score_frame.arm != "parent"].groupby("peptide_id", sort=False):
        parent_id, arm = group.parent_id.iloc[0], group.arm.iloc[0]
        parent_group = score_frame[score_frame.peptide_id == parent_id].set_index("model_number")
        group = group.set_index("model_number").sort_index()
        parent_group = parent_group.loc[group.index]
        passing = group[group.structure_pass.astype(bool)]
        common_pass = group.structure_pass.astype(bool) & parent_group.structure_pass.astype(bool)
        qualified_group, qualified_parent = group[common_pass], parent_group[common_pass]
        gain = qualified_parent.binding_score - qualified_group.binding_score
        pass_sequence = len(passing) >= 2
        med_gain = float(qualified_parent.binding_score.median() - qualified_group.binding_score.median()) if len(qualified_group) >= 2 else np.nan
        improved = bool(len(qualified_group) >= 2 and med_gain > .05 and int((gain > 0).sum()) >= 2)
        summaries.append({"peptide_id": ident, "parent_id": parent_id, "arm": arm,
                          "replicate": int(group.replicate.iloc[0]), "sequence": group.sequence.iloc[0],
                          "mask_positions": group.mask_positions.iloc[0], "mutations": group.mutations.iloc[0],
                          "anchor_positions_mutated": ";".join(str(pos) for pos in [3, 7, 10] if pos in [int(p) for p in str(group.mask_positions.iloc[0]).split(";")]),
                          "ppl_score": float(group.ppl_score.iloc[0]), "quality_pass_models": len(passing),
                          "common_passing_parent_child_models": len(qualified_group),
                          "sequence_pass": pass_sequence, "graphpep_median": float(group.binding_score.median()),
                          "predicted_parent_graphpep_median": float(parent_group.binding_score.median()),
                          "median_gain_vs_predicted_parent": med_gain,
                          "unfiltered_median_gain_vs_parent": float(parent_group.binding_score.median() - group.binding_score.median()),
                          "matched_models_improved": int((gain > 0).sum()), "predefined_improvement": improved,
                          "all_3_models_improved": bool(improved and int((gain > 0).sum()) == 3),
                          "iptm_mean": float(group.iptm.mean()), "peptide_plddt_mean": float(group.peptide_plddt.mean()),
                          "interface_pae_median": float(group.interface_pae.median()),
                          "contact_model_jaccard_mean": float(group.contact_model_jaccard.mean()),
                          "native_contact_recall_mean": float(group.native_contact_recall.mean()),
                          "peptide_CA_RMSD_mean_A": float(group.peptide_CA_RMSD_after_receptor_alignment_A.mean()),
                          "clash_models": int((group.clashing_atom_pairs_lt_1p5A > 0).sum())})
    seq = pd.DataFrame(summaries)
    seq.to_csv(OUT / "analysis/sequence_summary.csv", index=False)
    arm_rows = []
    for (parent_id, arm), group in seq.groupby(["parent_id", "arm"], sort=False):
        values = score_frame[(score_frame.parent_id == parent_id) & (score_frame.arm == arm)]
        arm_rows.append({"parent_id": parent_id, "arm": arm, "n_children": len(group),
                         "quality_comparable_children": int(group.median_gain_vs_predicted_parent.notna().sum()),
                         "structure_pass_sequences": int(group.sequence_pass.sum()),
                         "structure_pass_models": int(values.structure_pass.sum()), "n_models": len(values),
                         "children_improved": int(group.predefined_improvement.sum()),
                         "children_improved_all_3_models": int(group.all_3_models_improved.sum()),
                         "mean_child_median_score": float(group.graphpep_median.mean()),
                         "mean_child_median_gain_vs_parent": float(group.median_gain_vs_predicted_parent.mean())})
    arms = pd.DataFrame(arm_rows)
    arms.to_csv(OUT / "analysis/arm_summary.csv", index=False)
    native_scores = pd.read_csv(OUT / "native_scores.csv")
    parent_stats = []
    for parent in metadata["parents"]:
        native = native_scores[native_scores.peptide_id == parent["peptide_id"]].iloc[0]
        predicted = score_frame[score_frame.peptide_id == parent["peptide_id"]]
        parent_stats.append({"parent_id": parent["peptide_id"], "sequence": parent["sequence"],
                             "source_pdb": parent["pdb"], "native_coordinate_peptide_sequence": parent["native_peptide_coordinate_sequence"],
                             "native_unresolved_peptide_positions": parent["native_peptide_unresolved_positions"],
                             "native_graphpep_score": float(native.binding_score), "native_graphpep_contacts": int(native.n_contacts),
                             "predicted_parent_median_score": float(predicted.binding_score.median()),
                             "predicted_parent_structure_pass": int(predicted.structure_pass.sum()),
                             "predicted_parent_model_iptm_mean": float(predicted.iptm.mean()),
                             "predicted_parent_peptide_plddt_mean": float(predicted.peptide_plddt.mean())})
    write_json(OUT / "analysis/parent_summary.json", {"parents": parent_stats})
    source_hashes = {str(gp.GP / name): digest(gp.GP / name) for name in ["GraphPep.sh", "bin/core.py", "bin/data.py"]}
    if source_hashes != metadata["native_graphpep_source_hashes"]:
        raise ValueError("Native GraphPep source hash changed")
    if len(model_rows) != 117 or len(score_frame) != 117 or score_frame.loc[score_frame.graphpep_status == "ok", "binding_score"].isna().any():
        raise ValueError(f"Incomplete panel scores: structures={len(model_rows)}, scores={len(score_frame)}")
    verification = {"parent_count": 3, "sequence_count": 39, "predicted_structure_count": 117,
                    "graphpep_predicted_structure_scores": 117, "graphpep_native_parent_scores": 3,
                    "finite_graphpep_scores": int(score_frame.binding_score.notna().sum()),
                    "no_contact_structures": int((score_frame.graphpep_status == "no_contact").sum()),
                    "all_chain_sequences_verified": True, "all_confidence_dimensions_verified": True,
                    "all_graphpep_formulas_and_contacts_verified": True,
                    "all_within_arm_mutations_verified": True,
                    "graphpep_source_hashes_unchanged": True, "graphpep_source_hashes": source_hashes,
                    "new_random_seed": SEED}
    write_json(OUT / "analysis/verification.json", verification)
    write_report(metadata, seq, arms, parent_stats, verification)


def write_report(metadata: dict, seq: pd.DataFrame, arms: pd.DataFrame, parents: list[dict], verification: dict) -> None:
    def arm(parent_id: str, name: str) -> pd.Series:
        return arms[(arms.parent_id == parent_id) & (arms.arm == name)].iloc[0]
    lines = ["# MDM2多父代反馈位点对照与新种子复验", "",
             "**设计：** 复验1YCR父代并扩展到两个额外MDM2–肽PDB父代。每个父代比较实验结构GraphPep反馈、AFM父代结构GraphPep反馈及同预算随机mask；每臂3条、每条3个AFM参数模型。",
             "", "## 父代和结构", "", "|父代|来源|完整肽序列|实验结构解析肽|native GraphPep|预测父代结构通过|",
             "|---|---|---|---|---:|---:|"]
    for parent in parents:
        lines.append(f"|{parent['parent_id']}|{parent['source_pdb']}|`{parent['sequence']}`|`{parent['native_coordinate_peptide_sequence']}`|{parent['native_graphpep_score']:.4f}|{parent['predicted_parent_structure_pass']}/3|")
    lines += ["", "实验结构中的缺失肽末端残基使用完整PDB序列参与PepMLM和AFM，实验坐标及原生GraphPep反馈只使用实际解析残基；报告保留这一边界。新父代分别来自[RCSB 3EQS](https://www.rcsb.org/structure/3EQS)和[RCSB 3G03](https://www.rcsb.org/structure/3G03)。", "",
              "## 预先固定的比较标准", "",
              "新seed未松弛p53父代有两个模型各出现一处<1.5 Å的重原子重叠，因此在生成任何子代前统一加入ColabFold官方Amber有约束能量最小化。所有父代及子代都使用同样参数：200次以内/轮，最多3轮，stiffness=10，tolerance=2.39，CPU执行。原始结构和原始父代评分保留。结构置信度仍是AFM原始预测值，Amber不会重新估计ipTM或pLDDT。", "",
              "四个实验臂均替换2个位点、每个父代每臂生成3条、PepMLM top-k=3、无PPL筛选。随机臂从同一反馈来源记录到接触的位置均匀抽2位，与对应的native-guided或predicted-guided策略匹配。实验反馈来源为PDB原生复合物；预测反馈来自通过质量过滤的AFM父代模型（至少2/3），逐残基置信度按接触数量加权。AFM使用同一个新seed、相同受体MSA和3个参数模型。", "",
              "逐序列改善门槛：父代与子代至少有2个相同参数模型同时通过ipTM≥0.7、肽pLDDT≥70、界面PAE≤10 Å、最小重原子距离≥1.5 Å且存在界面接触；在这些通过过滤的配对模型上，中位GraphPep分数降低>0.05；并有至少2个配对模型方向一致。表中的GraphPep中位数保留全部3个模型用于诊断，增益只用通过质量过滤的配对模型计算，不足2个记为NA。GraphPep不是实验亲和力。", "",
              "GraphPep分数降低>0.05只是此小实验的操作性门槛，不对应经校准的亲和力变化或统计显著性。", "",
              "## 父代反馈一致性与选位", "", "|父代|native选位|预测反馈选位|共有反馈位点|Pearson|Spearman|", "|---|---|---|---:|---:|---:|"]
    for parent in metadata["parents"]:
        native_sites = ";".join(f"{parent['sequence'][pos-1]}{pos}" for pos in parent["native_guided_positions"])
        predicted_sites = ";".join(f"{parent['sequence'][pos-1]}{pos}" for pos in parent["predicted_guided_positions"])
        pearson = parent["native_predicted_feedback_pearson"]
        spearman = parent["native_predicted_feedback_spearman"]
        lines.append(f"|{parent['peptide_id']}|{native_sites}|{predicted_sites}|{parent['native_predicted_feedback_common_positions']}|{pearson if pearson is not None else 'NA'}|{spearman if spearman is not None else 'NA'}|")
    lines += ["", "## 每条子代", "", "|父代|选位策略|序列|PPL|GraphPep中位数|较预测父代增益|结构通过|达标|", "|---|---|---|---:|---:|---:|---:|---|"]
    for row in seq.to_dict("records"):
        gain = "NA" if pd.isna(row['median_gain_vs_predicted_parent']) else f"{row['median_gain_vs_predicted_parent']:+.4f}"
        lines.append(f"|{row['parent_id']}|{row['arm']}|`{row['sequence']}`|{row['ppl_score']:.2f}|{row['graphpep_median']:.4f}|{gain}|{row['quality_pass_models']}/3|{'是' if row['predefined_improvement'] else '否'}|")
    lines += ["", "## 每个父代分组", "", "|父代|策略|结构通过子代|达到改善标准|平均分数增益|", "|---|---|---:|---:|---:|"]
    for _, row in arms.iterrows():
        lines.append(f"|{row['parent_id']}|{row['arm']}|{row['structure_pass_sequences']}/3|{row['children_improved']}/3|{row['mean_child_median_gain_vs_parent']:+.4f}|")
    lines += ["", "## 结果解释", ""]
    for parent in PARENTS:
        ident = parent["peptide_id"]
        comparisons = []
        for guide, control in [("native_guided", "native_random"), ("predicted_guided", "predicted_random")]:
            g, r = arm(ident, guide), arm(ident, control)
            comparisons.append(f"{guide} {g.children_improved}/3，matched-random {r.children_improved}/3")
        lines.append(f"- **{ident}:** " + "；".join(comparisons) + ".")
    native_total = int(arms[arms.arm == "native_guided"].children_improved.sum())
    native_random = int(arms[arms.arm == "native_random"].children_improved.sum())
    predicted_total = int(arms[arms.arm == "predicted_guided"].children_improved.sum())
    predicted_random = int(arms[arms.arm == "predicted_random"].children_improved.sum())
    lines += ["", f"跨3个父代，native反馈指导 {native_total}/9 达到预设改善门槛，native对应随机组 {native_random}/9；AFM预测结构反馈指导 {predicted_total}/9，预测反馈对应随机组 {predicted_random}/9。",
              "一个父代仅3条子代，是低样本的初步对照，不做显著性声称；同一AFM权重下3个模型不是生物学重复。多个序列和模型上的一致趋势只构成计算证据。", "",
              "本试验评估同一MDM2结构域上的序列局部优化。不同父代的绝对GraphPep分数不能当作彼此亲和力排名；结构重现或GraphPep变好也不代表实验亲和力已提高。", "",
              "三条父代共享F3/W7/L10结合锚点。随机组可能更常改动锚点，需在`sequence_summary.csv`的`anchor_positions_mutated`中检查这一混杂因素；本轮还不能区分GraphPep选位收益和保留已知锚点的收益。逐残基反馈是模型接触置信度，不是物理结合能的逐残基分解。MDM2上的结果尚未迁移验证到NCAM1。", "",
              "## 核验与产物", "", f"完整性核验：{verification['parent_count']}个父代，{verification['sequence_count']}条序列，{verification['predicted_structure_count']}个AFM复合物，{verification['graphpep_predicted_structure_scores']}个预测结构GraphPep评分和{verification['graphpep_native_parent_scores']}个原生结构评分。链序列、置信度维数、GraphPep公式/接触数及GraphPep源码哈希均已独立检查。", "",
              "详情见`candidates.csv`、`analysis/model_results.csv`、`analysis/sequence_summary.csv`、`analysis/arm_summary.csv`、`analysis/parent_summary.json`和`analysis/verification.json`。"]
    (OUT / "report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def collect_native_contact_positions(parent: dict) -> set[tuple[int, int]]:
    rec = residue_atom_data(Path(parent["native_path"]), "A")
    pep = residue_atom_data(Path(parent["native_path"]), "B")
    pairs, _, _ = distances(rec, pep, 5.)
    return pairs  # native coordinate positions; metrics maps them to the AFM construct once


def score_offspring() -> None:
    metadata = json.loads((OUT / "parent_metadata.json").read_text())
    if metadata.get("stage") != "offspring_and_msas_ready":
        raise ValueError("Prepare offspring before analyzing")
    native_contacts = {p["peptide_id"]: collect_native_contact_positions(p) for p in metadata["parents"]}
    for parent in metadata["parents"]:
        parent["native_contact_positions"] = sorted([list(pair) for pair in native_contacts[parent["peptide_id"]]])
    write_json(OUT / "parent_metadata.json", metadata)
    analyze()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "score-parents", "generate", "analyze"])
    args = parser.parse_args()
    started = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage == "prepare":
        make_parent_inputs()
    elif args.stage == "score-parents":
        parent_structure_scores()
    elif args.stage == "generate":
        generate_children()
    elif args.stage == "analyze":
        score_offspring()
    print(f"Stage {args.stage} completed in {time.perf_counter()-started:.1f}s", flush=True)


if __name__ == "__main__":
    main()
