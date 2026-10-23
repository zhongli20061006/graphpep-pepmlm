#!/usr/bin/env python3
"""Disentangle peptide generation, predicted geometry and affinity-score validity."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import math
import os
from pathlib import Path
import re

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer

from affinity_reference import REFERENCES, SOURCE_DOI, SOURCE_URL
from affinity_reference import affinity_rank_correlation, paired_median_gain
import feedback_panel as panel

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/affinity_diagnostic_mdm2_20261006"
SEEDS = (20261011, 20261012)
TARGET = panel.TARGET
COMMON_REC = tuple(range(3, 85))  # canonical MDM2 27--108
COMMON_PEP = tuple(range(1, 12))
QUALITY_GATE = {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10.,
                "minimum_distance_min": 1.5}
NATIVE_SOURCES = {
    "PMI_WT": ROOT / "output/feedback_panel_mdm2_20261005/pdb/3EQS.cif",
    "PMI_N8A": OUT / "pdb/3LNZ.cif",
}
REF_BY_ID = {row.peptide_id: row for row in REFERENCES}
SCORE_OUTPUT_FIELDS = {
    "binding_score", "red_residue_ratio", "low_conf_contact_pair_ratio", "n_contacts",
    "min_conf", "mean_conf", "graphpep_status", "n_geometry_contact_pairs",
    "graphpep_workdir", "graphpep_seconds", "contact_feedback_path", "residue_feedback_path",
    "structure_sha256", "graphpep_source_signature", "contact_pairs_6A",
}


def digest(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def score_input_metadata(row: dict) -> dict:
    """A new coordinate input must never carry a previous structure's score."""
    return {key: value for key, value in row.items() if key not in SCORE_OUTPUT_FIELDS}


def configure_embedding_cache() -> None:
    """Keep this experiment independent of pre-existing unversioned embeddings."""
    cache = OUT / "graphpep_esm_cache"
    cache.mkdir(exist_ok=True)
    os.environ["GRAPHPEP_ESM_CACHE"] = str(cache)


def selected_atoms(residue) -> list:
    """Resolve alternative atoms by occupancy, blank/A preference, then altloc."""
    choices = {}
    for atom in residue.get_unpacked_list():
        if atom.element.upper() in {"H", "D"}:
            continue
        choices.setdefault(atom.name, []).append(atom)
    return [sorted(atoms, key=lambda atom: (-(atom.occupancy or 0.),
            0 if atom.altloc == " " else 1 if atom.altloc == "A" else 2, atom.altloc))[0]
            for atoms in choices.values()]


def write_complex(path: Path, receptor_residues: list, peptide_residues: list) -> None:
    """Write deterministic A/B heavy-atom PDBs without changing coordinates."""
    path.parent.mkdir(parents=True, exist_ok=True)
    serial, lines = 1, []
    for chain_id, residues in (("A", receptor_residues), ("B", peptide_residues)):
        for index, residue in enumerate(residues, 1):
            atoms = selected_atoms(residue)
            if not atoms:
                raise ValueError("Residue has no selected atoms")
            for atom in atoms:
                x, y, z = atom.coord
                element = atom.element.upper()
                name = atom.name.strip()
                atom_field = f" {name:<3s}" if len(name) < 4 and len(element) == 1 else f"{name:>4s}"
                occupancy = atom.occupancy if atom.occupancy is not None else 1.
                lines.append(f"ATOM  {serial:5d} {atom_field} {residue.resname:>3s} {chain_id}{index:4d}    "
                             f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy:6.2f}{atom.bfactor:6.2f}"
                             f"          {element:>2s}\n")
                serial += 1
        lines.append(f"TER   {serial:5d}      {residues[-1].resname:>3s} {chain_id}{len(residues):4d}\n")
        serial += 1
    path.write_text("".join(lines) + "END\n")


def observed_contact_count(receptor_residues: list, peptide_residues: list) -> int:
    def arrays(residues):
        xyz, indices = [], []
        for i, residue in enumerate(residues):
            for atom in selected_atoms(residue):
                xyz.append(atom.coord)
                indices.append(i)
        return np.asarray(xyz), np.asarray(indices)
    rec_xyz, rec_indices = arrays(receptor_residues)
    pep_xyz, pep_indices = arrays(peptide_residues[:11])
    distance = np.linalg.norm(rec_xyz[:, None] - pep_xyz[None, :], axis=2)
    i, j = np.where(distance <= 6.)
    return len(set(zip(rec_indices[i].tolist(), pep_indices[j].tolist())))


def prepare_native() -> list[dict]:
    rows = []
    for peptide_id, source in NATIVE_SOURCES.items():
        entry = source.stem.upper()
        expected = REF_BY_ID[peptide_id].sequence
        chemical_sequences = [re.sub(r"\s+", "", value)
                              for value in MMCIF2Dict(str(source))["_entity_poly.pdbx_seq_one_letter_code_can"]]
        if expected not in chemical_sequences or TARGET not in chemical_sequences:
            raise ValueError(f"Deposited chemical sequences differ from assay construct: {source}")
        model = MMCIFParser(QUIET=True).get_structure(entry, str(source))[0]
        receptors, peptides = {}, {}
        for chain in model:
            residues, sequence = panel.standard_residues(chain)
            if len(sequence) > 50:
                # Author numbering distinguishes unresolved terminal VV residues
                # that a sequence-only alignment can place at either 108 or 109.
                mapping = {i: residue.id[1] - 24 for i, residue in enumerate(residues, 1)}
                correct_sequence = all(1 <= pos <= len(TARGET) and
                                       TARGET[pos - 1] == panel.gp.THREE_TO_ONE[residues[i - 1].resname]
                                       for i, pos in mapping.items())
                if not correct_sequence or not set(COMMON_REC).issubset(mapping.values()):
                    raise ValueError(f"Native receptor does not contain common range: {entry}:{chain.id}")
                receptors[chain.id] = (residues, sequence, mapping)
            elif sequence and expected.startswith(sequence):
                if len(sequence) not in {11, 12}:
                    raise ValueError("Unexpected resolved peptide length")
                peptides[chain.id] = (residues, sequence)
        expected_copies = 1 if entry == "3EQS" else 8
        if len(receptors) != expected_copies or len(peptides) != expected_copies:
            raise ValueError(f"Unexpected crystallographic copy count in {entry}")
        for peptide_chain, (pep_res, pep_seq) in peptides.items():
            contacts = {chain: observed_contact_count(values[0], pep_res) for chain, values in receptors.items()}
            receptor_chain = max(contacts, key=contacts.get)
            # Validate against deposited alternating receptor/peptide chain order.
            if receptor_chain != chr(ord(peptide_chain) - 1):
                raise ValueError(f"Geometry pairing disagrees with deposited chain order: {contacts}")
            rec_res, rec_seq, mapping = receptors[receptor_chain]
            common_rec_res = [residue for i, residue in enumerate(rec_res, 1) if mapping[i] in COMMON_REC]
            for scope in ("resolved", "common"):
                receptor = rec_res if scope == "resolved" else common_rec_res
                peptide = pep_res if scope == "resolved" else pep_res[:11]
                sequence = pep_seq if scope == "resolved" else expected[:11]
                rec_sequence = rec_seq if scope == "resolved" else TARGET[2:84]
                ident = f"{entry}_{receptor_chain}{peptide_chain}_{scope}"
                path = OUT / "native_complexes" / f"{ident}.pdb"
                write_complex(path, receptor, peptide)
                rows.append({"peptide_id": peptide_id, "peptide_sequence": sequence,
                             "receptor_sequence": rec_sequence, "model_id": ident,
                             "complex_path": str(path), "source": "crystal", "scope": scope,
                             "pdb_id": entry, "original_receptor_chain": receptor_chain,
                             "original_peptide_chain": peptide_chain,
                             "chemical_sequence": expected, "kd_nM": REF_BY_ID[peptide_id].kd_nM,
                             "resolved_peptide_sequence": pep_seq,
                             "unresolved_peptide_positions": ";".join(map(str, range(len(pep_seq) + 1, 13))),
                             "resolved_receptor_canonical_positions": ";".join(map(str, mapping.values())),
                             "pairing_contact_count": contacts[receptor_chain]})
    return rows


def prepare() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "protocol.json").exists():
        verify_frozen()
        print("Existing protocol and inputs verified", flush=True)
        return
    rows = []
    for reference in REFERENCES:
        row = asdict(reference)
        if len(reference.sequence) != 12:
            raise ValueError("Assayed peptide must remain a 12-mer")
        parent = REF_BY_ID[reference.parent_id].sequence
        changes = [i for i, (a, b) in enumerate(zip(parent, reference.sequence), 1) if a != b]
        if changes != ([reference.position] if reference.position else []):
            raise ValueError("Sequence does not match documented mutation")
        row.update(kd_M=reference.kd_nM * 1e-9, sd_M=reference.sd_nM * 1e-9,
                   log10_kd_nM=math.log10(reference.kd_nM), source_url=SOURCE_URL,
                   source_doi=SOURCE_DOI, assay="SPR competition", receptor_construct="MDM2 25-109",
                   n_terminal="free amine", c_terminal="free acid", peptide_sequence=reference.sequence)
        rows.append(row)
    reference_path = OUT / "reference_affinity.csv"
    pd.DataFrame(rows).to_csv(reference_path, index=False)
    native_rows = prepare_native()
    native_manifest = OUT / "native_manifest.csv"
    pd.DataFrame(native_rows).to_csv(native_manifest, index=False)
    source_msa = panel.PRIOR / "results/POSCTRL_MDM2_P53.a3m"
    text = source_msa.read_text()
    msas, msa_stats = [], []
    for reference in REFERENCES:
        msa, stats = panel.make_shared_msa(text, TARGET, "ETFSDLWKLLPEN", reference.sequence)
        path = OUT / "afm_msa" / f"{reference.peptide_id}.a3m"
        path.parent.mkdir(exist_ok=True)
        path.write_text(msa)
        msas.append(path)
        msa_stats.append({"peptide_id": reference.peptide_id, **stats})
    if len({row["receptor_alignment_sha256"] for row in msa_stats}) != 1:
        raise ValueError("Unequal receptor alignments")
    (OUT / "analysis").mkdir(exist_ok=True)
    sources = [Path(__file__), ROOT / "pipeline/affinity_reference.py", ROOT / "pipeline/feedback_panel.py",
               ROOT / "pipeline/step3_graphpep.py", ROOT / "pipeline/shared_receptor_msa.py",
               ROOT / "pipeline/relax_feedback_panel.py", panel.gp.ROOT / "pipeline/graphpep_patches/resmap_compat.py",
               ROOT / "scripts/run_affinity_diagnostic.sh", reference_path, native_manifest, source_msa]
    sources += list(NATIVE_SOURCES.values()) + msas
    sources += [Path(row["complex_path"]) for row in native_rows]
    sources += [panel.gp.GP / "GraphPep.sh"]
    sources += sorted((panel.gp.GP / "bin").rglob("*.py"))
    sources += sorted((panel.gp.GP / "bin").glob("*.awk"))
    sources += [panel.gp.GP / "bin/itpot.csv"]
    sources += [panel.gp.GP / "ckpt/model_param.ckpt",
                panel.gp.GP / "models/esm2_t33_650M_UR50D.pt",
                panel.gp.GP / "models/esm2_t33_650M_UR50D-contact-regression.pt"]
    protocol = {
        "state": "frozen_before_any_benchmark_score", "source": SOURCE_URL,
        "receptor_sequence": TARGET, "seeds": SEEDS, "models": [1, 2, 3], "templates": False,
        "recycles": 3, "early_stop_tolerance": 0., "max_seq": 256, "max_extra_seq": 512,
        "relaxation": panel.COMPARISON_PROTOCOL["relaxation"], "quality_gate": QUALITY_GATE,
        "minimum_passing_models_per_seed": 2,
        "primary_correlation": "Spearman(GraphPep score, log10(Kd_nM)); positive is correct",
        "pair_gain": "median of matched parent-minus-variant model score differences",
        "strong_direction_cases": [row.peptide_id for row in REFERENCES if row.position and
                                   max(row.kd_nM / REF_BY_ID[row.parent_id].kd_nM,
                                       REF_BY_ID[row.parent_id].kd_nM / row.kd_nM) >= 3],
        "continuation_rule": "Each seed: overall rho >=0.6; both scaffold rhos >0; >=6/7 clear directions correct, including both stronger variants; all seven comparable",
        "native_common_receptor_canonical_positions": COMMON_REC,
        "native_common_receptor_pdb_range": "27-108",
        "native_common_peptide_positions": COMMON_PEP,
        "native_copies_retained": {"3EQS": 1, "3LNZ": 8},
        "alternate_atoms": "highest occupancy; blank then A then lexicographic altloc",
        "graphpep_embedding_cache": str(OUT / "graphpep_esm_cache"),
        "exclude_1YCR_as_exact_affinity_control": "different deposited peptide/receptor chemical constructs",
        "assay_conditions": {"temperature": "room temperature; no exact numeric value provided",
                             "buffer": "10 mM HEPES, 150 mM NaCl, 0.005% P20, pH 7.4",
                             "receptor_nM": 50, "replicates": "at least three independent measurements",
                             "immobilized_peptide": "p53(15-29)"},
        "msa_stats": msa_stats,
        "source_hashes": {str(path): digest(path) for path in sources},
    }
    write_json(OUT / "protocol.json", protocol)
    print(f"Frozen {len(rows)} assayed sequences, {len(native_rows)//2} crystal complexes, {len(SEEDS)*3*len(rows)} AFM models", flush=True)


def verify_frozen() -> dict:
    protocol = json.loads((OUT / "protocol.json").read_text())
    for raw, expected in protocol["source_hashes"].items():
        if digest(Path(raw)) != expected:
            raise ValueError(f"Frozen input/source changed: {raw}")
    return protocol


def quality_metrics(reference: dict, seed: int, path: Path) -> dict:
    rec = panel.residue_atom_data(path, "A")
    pep = panel.residue_atom_data(path, "B")
    if rec[0] != TARGET or pep[0] != reference["sequence"]:
        raise ValueError(f"Predicted sequence mismatch: {path}")
    confidence = path.with_name(path.name.replace("_relaxed_", "_scores_")).with_suffix(".json")
    values = json.loads(confidence.read_text())
    pae, plddt = np.asarray(values["pae"]), np.asarray(values["plddt"])
    nr, np_ = len(TARGET), len(reference["sequence"])
    if pae.shape != (nr + np_, nr + np_) or plddt.shape != (nr + np_,):
        raise ValueError("Confidence shape differs from chain lengths")
    contacts, minimum, clashes = panel.distances(rec, pep, 6.)
    contact_rec = sorted({r - 1 for r, _ in contacts})
    pep_indices = list(range(nr, nr + np_))
    interface_pae = float(np.median(np.concatenate([
        pae[np.ix_(contact_rec, pep_indices)].ravel(),
        pae[np.ix_(pep_indices, contact_rec)].ravel()]))) if contacts else math.nan
    peptide_plddt = float(plddt[nr:].mean())
    failures = []
    if not math.isfinite(float(values["iptm"])) or float(values["iptm"]) < QUALITY_GATE["iptm_min"]:
        failures.append("low_or_nonfinite_iptm")
    if not math.isfinite(peptide_plddt) or peptide_plddt < QUALITY_GATE["peptide_plddt_min"]:
        failures.append("low_or_nonfinite_peptide_plddt")
    if not math.isfinite(interface_pae) or interface_pae > QUALITY_GATE["interface_pae_max"]:
        failures.append("high_or_missing_interface_pae")
    if minimum < QUALITY_GATE["minimum_distance_min"]: failures.append("heavyatom_clash")
    if not contacts: failures.append("no_interface_contacts")
    model_number = int(re.search(r"_model_(\d+)_seed_", path.name).group(1))
    return {**reference, "seed": seed, "model_number": model_number,
            "model_id": f"seed{seed}_model{model_number}",
            "complex_path": str(path), "receptor_sequence": TARGET, "scope": "full", "source": "AFM",
            "confidence_path": str(confidence), "iptm": float(values["iptm"]), "ptm": float(values["ptm"]),
            "peptide_plddt": peptide_plddt, "interface_pae": interface_pae,
            "minimum_heavyatom_distance": minimum, "clashing_atom_pairs_lt_1p5A": clashes,
            "contact_pairs_6A": len(contacts), "structure_pass": not failures,
            "quality_fail_reasons": ";".join(failures)}


def score_afm(seed: int) -> None:
    verify_frozen()
    configure_embedding_cache()
    references = pd.read_csv(OUT / "reference_affinity.csv", keep_default_na=False)
    rows = []
    directory = OUT / f"seed_{seed}" / "relaxed"
    for reference in references.to_dict("records"):
        paths = sorted(directory.glob(f"{reference['peptide_id']}_relaxed_rank_*_model_*_seed_{seed:03d}.pdb"))
        if len(paths) != 3:
            raise ValueError(f"Expected 3 models for {reference['peptide_id']} seed {seed}, found {len(paths)}")
        rows.extend(quality_metrics(reference, seed, path) for path in paths)
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "analysis" / f"model_metrics_seed{seed}.csv", index=False)
    scores = panel.run_graphpep(metrics, OUT / "analysis" / f"graphpep_seed{seed}")
    for row in scores.to_dict("records"):
        panel.score_bounds(row)
    scores.to_csv(OUT / "analysis" / f"model_results_seed{seed}.csv", index=False)


def score_native() -> None:
    verify_frozen()
    configure_embedding_cache()
    native_manifest = pd.read_csv(OUT / "native_manifest.csv", keep_default_na=False)
    # Shared coordinate inputs can occur for an already-common resolved copy.
    # The cache keys keep the scopes explicit, and no experimental atoms move.
    native_scores = panel.run_graphpep(native_manifest, OUT / "analysis/native_graphpep")
    native_scores["contact_pairs_6A"] = native_scores.n_geometry_contact_pairs
    for row in native_scores.to_dict("records"):
        row["contact_pairs_6A"] = row["n_geometry_contact_pairs"]
        panel.score_bounds(row)
    native_scores.to_csv(OUT / "analysis/native_results.csv", index=False)
    parser = PDBParser(QUIET=True)
    common_rows = []
    for seed in SEEDS:
        model_rows = pd.read_csv(OUT / "analysis" / f"model_results_seed{seed}.csv")
        for row in model_rows[model_rows.peptide_id.isin(NATIVE_SOURCES)].to_dict("records"):
            model = parser.get_structure("afm", row["complex_path"])[0]
            receptor, _ = panel.standard_residues(model["A"])
            peptide, _ = panel.standard_residues(model["B"])
            path = OUT / "predicted_common" / f"{row['peptide_id']}_{row['model_id']}_common.pdb"
            write_complex(path, [receptor[i - 1] for i in COMMON_REC], [peptide[i - 1] for i in COMMON_PEP])
            common_rows.append({**score_input_metadata(row), "complex_path": str(path), "scope": "common",
                                "full_complex_path": row["complex_path"],
                                "receptor_sequence": TARGET[2:84], "peptide_sequence": row["sequence"][:11]})
    common_scores = panel.run_graphpep(pd.DataFrame(common_rows), OUT / "analysis/predicted_common_graphpep")
    common_scores["contact_pairs_6A"] = common_scores.n_geometry_contact_pairs
    for row in common_scores.to_dict("records"):
        row["contact_pairs_6A"] = row["n_geometry_contact_pairs"]
        panel.score_bounds(row)
    common_scores.to_csv(OUT / "analysis/predicted_common_results.csv", index=False)


def sequence_summaries(models: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (peptide_id, seed), group in models.groupby(["peptide_id", "seed"]):
        passing = group[group.structure_pass & group.binding_score.notna()]
        reference = REF_BY_ID[peptide_id]
        rows.append({**asdict(reference), "seed": int(seed), "n_models": len(group),
                     "passing_models": len(passing), "median_score_all": float(group.binding_score.median()),
                     "median_score_passing": float(passing.binding_score.median()) if len(passing) >= 2 else math.nan,
                     "iptm_min": float(group.iptm.min()), "peptide_plddt_min": float(group.peptide_plddt.min()),
                     "graphpep_score_range": float(group.binding_score.max() - group.binding_score.min())})
    return pd.DataFrame(rows)


def direction_summary(models: pd.DataFrame, protocol: dict) -> pd.DataFrame:
    rows = []
    for seed in SEEDS:
        subset = models[models.seed == seed]
        for reference in REFERENCES:
            if not reference.position:
                continue
            parent = subset[subset.peptide_id == reference.parent_id].set_index("model_number")
            variant = subset[subset.peptide_id == reference.peptide_id].set_index("model_number")
            common = [i for i in (1, 2, 3) if bool(parent.loc[i, "structure_pass"]) and
                      bool(variant.loc[i, "structure_pass"]) and
                      math.isfinite(float(parent.loc[i, "binding_score"])) and
                      math.isfinite(float(variant.loc[i, "binding_score"]))]
            expected_gain_sign = 1 if reference.kd_nM < REF_BY_ID[reference.parent_id].kd_nM else -1
            gain = paired_median_gain(parent.loc[common, "binding_score"].tolist(),
                                     variant.loc[common, "binding_score"].tolist()) if len(common) >= 2 else math.nan
            model_gains = (parent.loc[common, "binding_score"] - variant.loc[common, "binding_score"]).tolist()
            rows.append({"peptide_id": reference.peptide_id, "parent_id": reference.parent_id, "seed": seed,
                         "mutation": reference.mutation, "kd_ratio_variant_over_parent": reference.kd_nM / REF_BY_ID[reference.parent_id].kd_nM,
                         "expected_direction": "stronger" if expected_gain_sign == 1 else "weaker",
                         "primary_direction_case": reference.peptide_id in protocol["strong_direction_cases"],
                         "common_passing_models": len(common), "paired_median_gain": gain,
                         "correct_direction": math.isfinite(gain) and gain * expected_gain_sign > 0,
                         "models_correct_direction": sum(value * expected_gain_sign > 0 for value in model_gains)})
    return pd.DataFrame(rows)


def correlations(summaries: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for seed in SEEDS:
        for scope in ("all", "PMI", "p53"):
            subset = summaries[summaries.seed == seed]
            if scope != "all":
                subset = subset[subset.family == scope]
            for metric in ("median_score_all", "median_score_passing"):
                valid = subset[subset[metric].notna()]
                rho = affinity_rank_correlation(valid.kd_nM.tolist(), valid[metric].tolist()) if len(valid) >= 3 else math.nan
                rows.append({"seed": seed, "family": scope, "score_summary": metric,
                             "n_comparable_peptides": len(valid), "spearman_vs_log10_kd": rho})
    return pd.DataFrame(rows)


def native_comparison() -> pd.DataFrame:
    native = pd.read_csv(OUT / "analysis/native_results.csv")
    predicted = pd.read_csv(OUT / "analysis/predicted_common_results.csv")
    rows = []
    for label, subset in (("crystal_common", native[native.scope == "common"]),
                          ("crystal_resolved", native[native.scope == "resolved"]),
                          ("AFM_common", predicted)):
        groupings = ["peptide_id", "seed"] if label == "AFM_common" else ["peptide_id"]
        for keys, group in subset.groupby(groupings):
            ident = keys[0] if isinstance(keys, tuple) else keys
            seed = int(keys[1]) if label == "AFM_common" else 0
            rows.append({"peptide_id": ident, "source_scope": label, "seed": seed,
                         "kd_nM_full_chemical_12mer": REF_BY_ID[ident].kd_nM, "n_structural_copies": len(group),
                         "median_graphpep": float(group.binding_score.median()),
                         "min_graphpep": float(group.binding_score.min()), "max_graphpep": float(group.binding_score.max()),
                         "correct_affinity_unit": "full chemical peptide; NOT affinity of cropped coordinate sequence"})
    comparisons = pd.DataFrame(rows)
    comparisons.to_csv(OUT / "analysis/native_vs_predicted.csv", index=False)
    rmsds = []
    for crystal in native[native.scope == "common"].to_dict("records"):
        native_rec = panel.residue_atom_data(Path(crystal["complex_path"]), "A")
        native_pep = panel.residue_atom_data(Path(crystal["complex_path"]), "B")
        for prediction in predicted[predicted.peptide_id == crystal["peptide_id"]].to_dict("records"):
            pred_rec = panel.residue_atom_data(Path(prediction["complex_path"]), "A")
            pred_pep = panel.residue_atom_data(Path(prediction["complex_path"]), "B")
            fixed = np.asarray([native_rec[4][i] for i in range(1, 83)])
            moving = np.asarray([pred_rec[4][i] for i in range(1, 83)])
            sup = SVDSuperimposer()
            sup.set(fixed, moving)
            sup.run()
            rotation, translation = sup.get_rotran()
            peptide_difference = np.asarray([pred_pep[4][i] @ rotation + translation - native_pep[4][i] for i in range(1, 12)])
            rmsds.append({"peptide_id": crystal["peptide_id"], "native_model_id": crystal["model_id"],
                          "seed": prediction["seed"], "model_number": prediction["model_number"],
                          "receptor_CA_RMSD_A": float(sup.get_rms()),
                          "peptide_CA_RMSD_after_receptor_alignment_A": float(np.sqrt(np.mean(np.sum(peptide_difference**2, axis=1))))})
    pd.DataFrame(rmsds).to_csv(OUT / "analysis/native_afm_rmsd.csv", index=False)
    return comparisons


def feedback_diagnostic(models: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for seed in SEEDS:
        for parent_id in ("PMI_WT", "P53_WT"):
            parent = REF_BY_ID[parent_id]
            frames = []
            for row in models[(models.seed == seed) & (models.peptide_id == parent_id) & models.structure_pass].to_dict("records"):
                if row["graphpep_status"] == "ok":
                    frames.append(pd.read_csv(row["contact_feedback_path"]))
            contacts = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
            profile = contacts.groupby("peptide_position").CONFIDENT_SCORES.agg(["count", "mean"]) if len(contacts) else pd.DataFrame()
            eligible = [(i, float(profile.loc[i, "mean"])) for i in profile.index if i not in {3, 7, 10}]
            selected = {i for i, _ in sorted(eligible, key=lambda value: (value[1], value[0]))[:2]}
            for reference in REFERENCES:
                if reference.parent_id != parent_id or not reference.position:
                    continue
                position = reference.position
                present = position in profile.index
                rows.append({"parent_id": parent_id, "seed": seed, "peptide_id": reference.peptide_id,
                             "mutation": reference.mutation, "position_in_peptide": position,
                             "residue": parent.sequence[position - 1], "contact_observed": present,
                             "contact_count_all_passing_models": int(profile.loc[position, "count"]) if present else 0,
                             "mean_contact_confidence": float(profile.loc[position, "mean"]) if present else math.nan,
                             "protected_anchor": position in {3, 7, 10},
                             "selected_by_previous_low_confidence_rule": position in selected,
                             "kd_ratio_variant_over_parent": reference.kd_nM / parent.kd_nM,
                             "experimental_direction": "stronger" if reference.kd_nM < parent.kd_nM else "weaker"})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "analysis/feedback_vs_mutation.csv", index=False)
    return result


def plot(summaries: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharex=True, sharey=True)
    for axis, seed in zip(axes, SEEDS):
        subset = summaries[summaries.seed == seed]
        for family, color in (("PMI", "#2376ab"), ("p53", "#cc6e29")):
            group = subset[subset.family == family]
            axis.scatter(np.log10(group.kd_nM), group.median_score_all, color=color, label=family)
            for row in group.to_dict("records"):
                axis.annotate(row["mutation"], (math.log10(row["kd_nM"]), row["median_score_all"]),
                              xytext=(3, 4), textcoords="offset points", fontsize=7)
        axis.set_title(f"AFM seed {seed}; all-model medians")
        axis.set_xlabel("log10 experimental Kd (nM); lower = stronger")
        axis.grid(alpha=.2)
        axis.legend()
    axes[0].set_ylabel("GraphPep score; lower = better pose score")
    fig.tight_layout()
    fig.savefig(OUT / "analysis/affinity_scatter.png", dpi=180)
    plt.close(fig)


def analyze() -> None:
    protocol = verify_frozen()
    models = pd.concat([pd.read_csv(OUT / "analysis" / f"model_results_seed{seed}.csv") for seed in SEEDS], ignore_index=True)
    if len(models) != 60 or models.duplicated(["peptide_id", "seed", "model_number"]).any():
        raise ValueError("Expected sixty unique sequence/seed/model tuples")
    for row in models.to_dict("records"):
        panel.score_bounds(row)
    summaries = sequence_summaries(models)
    directions = direction_summary(models, protocol)
    corr = correlations(summaries)
    native = native_comparison()
    feedback = feedback_diagnostic(models)
    models.to_csv(OUT / "analysis/model_results.csv", index=False)
    summaries.to_csv(OUT / "analysis/sequence_summary.csv", index=False)
    directions.to_csv(OUT / "analysis/mutation_directions.csv", index=False)
    corr.to_csv(OUT / "analysis/rank_correlations.csv", index=False)
    checks = []
    for seed in SEEDS:
        current = corr[(corr.seed == seed) & (corr.score_summary == "median_score_passing")].set_index("family")
        cases = directions[(directions.seed == seed) & directions.primary_direction_case]
        stronger = cases[cases.expected_direction == "stronger"]
        ready = bool(current.loc["all", "n_comparable_peptides"] == 10 and
                     current.loc["all", "spearman_vs_log10_kd"] >= .6 and
                     current.loc["PMI", "spearman_vs_log10_kd"] > 0 and
                     current.loc["p53", "spearman_vs_log10_kd"] > 0 and
                     len(cases) == 7 and bool((cases.common_passing_models >= 2).all()) and
                     int(cases.correct_direction.sum()) >= 6 and bool(stronger.correct_direction.all()))
        checks.append({"seed": seed, "continuation_rule_pass": ready,
                       "clear_directions_correct": int(cases.correct_direction.sum()),
                       "clear_directions_comparable": int((cases.common_passing_models >= 2).sum())})
    decision = bool(all(row["continuation_rule_pass"] for row in checks))
    verify_frozen()
    verification = {"unique_reference_peptides": len(REFERENCES), "afm_models": len(models),
                    "passing_afm_models": int(models.structure_pass.sum()),
                    "finite_graphpep_scores": int(models.binding_score.notna().sum()),
                    "crystal_complexes": int((pd.read_csv(OUT / "analysis/native_results.csv").scope == "common").sum()),
                    "protocol_sha256": digest(OUT / "protocol.json"), "frozen_sources_unchanged": True,
                    "seed_checks": checks, "continue_current_affinity_objective": decision,
                    "new_peptides_generated": 0, "experimental_Kd_of_new_children_measured": False}
    write_json(OUT / "analysis/verification.json", verification)
    plot(summaries)
    write_report(summaries, directions, corr, native, feedback, verification)
    print(json.dumps(verification, ensure_ascii=False, indent=2), flush=True)


def write_report(summaries: pd.DataFrame, directions: pd.DataFrame, corr: pd.DataFrame,
                 native: pd.DataFrame, feedback: pd.DataFrame, verification: dict) -> None:
    lines = ["# MDM2：实测亲和力、结构来源与反馈选位诊断", "",
             "本轮绕过 PepMLM，使用十条已知的 12-aa L 肽，不生成新子代。Kd 来自同一篇论文的同一 SPR 竞争实验与 synthetic MDM2(25–109) 构建体；没有把 MDMX、不同肽长度或其他 assay 的数值混入。",
             f"数据来源：[Li et al. 2010 Table 1/3 与方法]({SOURCE_URL})。肽为游离 N/C 端，Kd/SD 由 M 换算为 nM。", "",
             "## 已知亲和力与 AFM→GraphPep", "",
             "|肽|序列|实测 Kd±SD (nM)|seed 20261011 全模型中位分数|seed 20261012 全模型中位分数|通过模型数（两个 seed）|",
             "|---|---|---:|---:|---:|---:|"]
    for reference in REFERENCES:
        values = summaries[summaries.peptide_id == reference.peptide_id].set_index("seed")
        lines.append(f"|{reference.peptide_id}|`{reference.sequence}`|{reference.kd_nM:g} ± {reference.sd_nM:g}|"
                     f"{values.loc[SEEDS[0], 'median_score_all']:.4f}|{values.loc[SEEDS[1], 'median_score_all']:.4f}|"
                     f"{int(values.passing_models.sum())}/6|")
    lines += ["", "GraphPep 越低与 Kd 越低方向一致，因此对 log10(Kd) 的相关系数应为正。两个 seed、六个 AFM 模型是计算敏感性检查，不是生物学重复。全模型结果与质量过滤结果分别报告，未通过结构门槛的弱结合肽不能被静默丢掉。", "",
              "|seed|范围|聚合口径|可比较肽数|Spearman 与 log10(Kd)|", "|---|---|---|---:|---:|"]
    for row in corr.to_dict("records"):
        lines.append(f"|{row['seed']}|{row['family']}|{row['score_summary']}|{row['n_comparable_peptides']}|{row['spearman_vs_log10_kd']:.3f}|")
    lines += ["", "## 单点突变方向", "",
              "只将实测变化至少三倍的七项作为主要方向检查；K24A 的变化不到三倍，仅作描述。增益固定为逐模型配对的父代分数减突变分数之差的中位数，正数表示 GraphPep 改善。", "",
              "|seed|突变|实测方向|Kd 比值(突变/父代)|配对中位增益|方向正确|可比模型|", "|---|---|---|---:|---:|---|---:|"]
    for row in directions.to_dict("records"):
        lines.append(f"|{row['seed']}|{row['peptide_id']}|{row['expected_direction']}|{row['kd_ratio_variant_over_parent']:.3g}|"
                     f"{row['paired_median_gain']:+.4f}|{'是' if row['correct_direction'] else '否/不可比较'}|{row['common_passing_models']}/3|")
    lines += ["", "## 实验结构与预测结构", "",
              "只比较化学序列匹配的 PMI(3EQS) 与 N8A-PMI(3LNZ)，保留 3LNZ 全部八个晶体副本。所有副本共同坐标为 MDM2 27–108 和肽 1–11；native 原始重原子坐标不移动，AFM 完整 85+12-aa 预测后按同一区域裁剪。原样解析范围评分另列，不能把裁剪序列当成其有 Kd 的截短肽。", "",
              "|序列|结构来源/范围|seed|结构副本数|中位分数|最小—最大分数|", "|---|---|---|---:|---:|---:|"]
    for row in native.to_dict("records"):
        lines.append(f"|{row['peptide_id']}|{row['source_scope']}|{row['seed']}|{row['n_structural_copies']}|"
                     f"{row['median_graphpep']:.4f}|{row['min_graphpep']:.4f} — {row['max_graphpep']:.4f}|")
    lines += ["", "[3EQS](https://www.rcsb.org/structure/3EQS) 和 [3LNZ](https://www.rcsb.org/structure/3LNZ) 的晶体副本不是独立亲和力测量，仅两个化学序列的对照不能证明普遍排名准确率。1YCR 的化学构建体与 p53(17–28) 不同，未给它的 native score 赋予 440 nM。原始坐标差异与 CA RMSD 见 `analysis/native_afm_rmsd.csv`。", "",
              "## 低置信度选位规则", "",
              "在保护 F3/W7/L10 的旧规则下，仅接触被观察到的非锚点参与选位；按父代通过质量门槛模型汇总接触置信度，取最低两位。以下列出本基准中实验变强的三项，K24A 仍为小幅变化。", "",
              "|seed|突变|有接触反馈|接触数|平均置信度|旧规则会选中|Kd 比值|", "|---|---|---|---:|---:|---|---:|"]
    for row in feedback[feedback.experimental_direction == "stronger"].to_dict("records"):
        lines.append(f"|{row['seed']}|{row['peptide_id']}|{row['contact_observed']}|{row['contact_count_all_passing_models']}|"
                     f"{row['mean_contact_confidence']:.3f}|{row['selected_by_previous_low_confidence_rule']}|{row['kd_ratio_variant_over_parent']:.3g}|")
    lines += ["", "## 冻结判据与下一步", "",
              "本轮评分前的操作性门槛：每个 seed 的总体 rho≥0.6、两个家族 rho 均>0；七项主要方向至少六项正确，含两个变强对照，全部七项可比较。此门槛不代表统计显著性或普遍适用性。", ""]
    for check in verification["seed_checks"]:
        lines.append(f"- seed {check['seed']}：主要方向正确 {check['clear_directions_correct']}/7，可比较 {check['clear_directions_comparable']}/7；整体门槛 {'通过' if check['continuation_rule_pass'] else '未通过'}。")
    if verification["continue_current_affinity_objective"]:
        lines += ["", "这套小基准支持继续检验本体系内的选位策略，但不证明新序列实际亲和力改善。下一轮应分别比较位点选择与氨基酸替换建议，并保持候选池/算力匹配。"]
    else:
        lines += ["", "当前 GraphPep 亲和力迭代目标未通过预定门槛，本轮停止沿该目标生成新子代。GraphPep 可保留为结构/接触评价工具；先建立独立亲和力目标及其正反基准，再继续生成闭环。该结果不能据此认定 PepMLM 有错，也不能仅凭高 AFM 置信度排除结构误差。"]
    lines += ["", "原论文指出 N8A 的收益涉及非接触位点的螺旋倾向。当前接触置信度不包含已验证的突变收益或游离肽构象熵分解；这需要单独的目标与实验检查。", "",
              "## 核验与产物", "",
              f"AFM 模型 {verification['afm_models']} 个，通过门槛 {verification['passing_afm_models']} 个；有效 GraphPep 分数 {verification['finite_graphpep_scores']} 个；实验结构副本 {verification['crystal_complexes']} 个。冻结数据/脚本/GraphPep 源码/共同 MSA 哈希已核对。", "",
              "原始逐模型数据 `analysis/model_results.csv`；双 seed 逐肽统计 `analysis/sequence_summary.csv`；方向检查 `analysis/mutation_directions.csv`；选位诊断 `analysis/feedback_vs_mutation.csv`；核验 `analysis/verification.json`。本轮不训练模型、不测量新子代 Kd，也没有计算 ΔG(binding)。", "",
              "![实测亲和力与 GraphPep](analysis/affinity_scatter.png)"]
    (OUT / "report_zh.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "score-afm", "score-native", "analyze"))
    parser.add_argument("--seed", type=int, choices=SEEDS)
    args = parser.parse_args()
    if args.stage == "prepare": prepare()
    elif args.stage == "score-afm":
        if args.seed is None: parser.error("score-afm requires --seed")
        score_afm(args.seed)
    elif args.stage == "score-native": score_native()
    else: analyze()


if __name__ == "__main__":
    main()
