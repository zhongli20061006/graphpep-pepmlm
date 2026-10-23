#!/usr/bin/env python3
"""Assemble paired pilot diagnostics without promoting failed structures to evidence."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser


PILOT = Path(__file__).resolve().parent
ROOT = PILOT.parents[1]
PARSER = PDBParser(QUIET=True)


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def ca_coords(path: Path, chain_id: str) -> np.ndarray:
    chain = PARSER.get_structure(path.stem, str(path))[0][chain_id]
    return np.asarray([res["CA"].coord for res in chain if res.id[0] == " "], dtype=float)


def contact_set(path: Path, cutoff: float = 6.0) -> set[tuple[int, int]]:
    model = PARSER.get_structure(path.stem, str(path))[0]
    atom_data = {}
    for chain_id in ("A", "B"):
        coords, residue_ix = [], []
        for i, residue in enumerate(r for r in model[chain_id] if r.id[0] == " "):
            for atom in residue:
                if atom.element.upper() not in {"H", "D"}:
                    coords.append(atom.coord.astype(float))
                    residue_ix.append(i)
        atom_data[chain_id] = (np.asarray(coords), np.asarray(residue_ix, dtype=int))
    a, ai = atom_data["A"]
    b, bi = atom_data["B"]
    pairs = np.argwhere(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2) <= cutoff)
    return {(int(ai[i]), int(bi[j])) for i, j in pairs}


def receptor_aligned_peptide_rmsd(ref: Path, mobile: Path) -> float:
    ref_a, mob_a = ca_coords(ref, "A"), ca_coords(mobile, "A")
    ref_b, mob_b = ca_coords(ref, "B"), ca_coords(mobile, "B")
    if ref_a.shape != mob_a.shape or ref_b.shape != mob_b.shape:
        raise ValueError("cross-seed CA counts differ")
    mob_center, ref_center = mob_a.mean(axis=0), ref_a.mean(axis=0)
    covariance = (mob_a - mob_center).T @ (ref_a - ref_center)
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    aligned_b = (mob_b - mob_center) @ rotation + ref_center
    return float(np.sqrt(np.mean(np.sum((aligned_b - ref_b) ** 2, axis=1))))


def sign(value: float) -> str:
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def main() -> int:
    config = json.loads((PILOT / "run_config_frozen_v1.json").read_text(encoding="utf-8"))
    gate = json.loads((PILOT / "structure_gate_audit_v1.json").read_text(encoding="utf-8"))
    gates = {(int(r["seed"]), r["complex_id"]): r for r in gate["rows"]}
    graphpep = {(int(r["model_id"].removeprefix("seed")), r["peptide_id"]): r for r in read_csv(PILOT / "scores/graphpep.csv")}
    foldx = {(int(r["seed"]), r["peptide_id"]): r for r in read_csv(PILOT / "scores/foldx.csv")}
    esm = {r["sequence_id"]: r for r in read_csv(PILOT / "scores/esm_pll.csv")}
    variants = ["BCLW_Bim_A2eT_I3dL_gain", "BCLW_Bim_A2eT_I2dM_loss"]
    parent_id = "BCLW_Bim_A2eT_parent"
    paired = []
    for seed in config["prediction"]["seeds"]:
        parent_g = float(graphpep[(seed, parent_id)]["binding_score"])
        parent_f = float(foldx[(seed, parent_id)]["interaction_energy_kcal_mol"])
        for variant in variants:
            g = float(graphpep[(seed, variant)]["binding_score"])
            f = float(foldx[(seed, variant)]["interaction_energy_kcal_mol"])
            parent_contacts = int(gates[(seed, parent_id)]["unique_contact_residue_pairs_le_6A"])
            variant_contacts = int(gates[(seed, variant)]["unique_contact_residue_pairs_le_6A"])
            paired.append({
                "seed": seed,
                "parent_id": parent_id,
                "variant_id": variant,
                "structure_quality_pass_for_both": bool(gates[(seed, parent_id)]["quality_pass"] and gates[(seed, variant)]["quality_pass"]),
                "GraphPep_parent_raw": parent_g,
                "GraphPep_variant_raw": g,
                "GraphPep_parent_minus_variant_proxy_gain": parent_g - g,
                "GraphPep_proxy_sign": sign(parent_g - g),
                "FoldX_parent_interaction_kcal_mol": parent_f,
                "FoldX_variant_interaction_kcal_mol": f,
                "FoldX_parent_minus_variant_proxy_gain_kcal_mol": parent_f - f,
                "FoldX_proxy_sign": sign(parent_f - f),
                "contact_pairs_parent": parent_contacts,
                "contact_pairs_variant": variant_contacts,
                "contact_pairs_variant_minus_parent": variant_contacts - parent_contacts,
                "pair_direction_evaluable_under_frozen_rule": False,
                "failure_reason": "both structures fail the frozen ipTM, peptide-pLDDT, and interface-PAE gates",
            })
    esm_parent = float(esm[parent_id]["sum_masked_log_likelihood"])
    esm_pairs = []
    for variant in variants:
        value = float(esm[variant]["sum_masked_log_likelihood"])
        esm_pairs.append({"parent_id": parent_id, "variant_id": variant, "parent_sum_PLL": esm_parent,
                          "variant_sum_PLL": value, "variant_minus_parent_sequence_prior": value - esm_parent,
                          "interpretation": "sequence plausibility only; not affinity or interface score"})

    pose = []
    for ident in [parent_id, *variants]:
        p08 = next((PILOT / "seed_20261008/relaxed").glob(f"{ident}_relaxed_rank_001_*.pdb"))
        p09 = next((PILOT / "seed_20261009/relaxed").glob(f"{ident}_relaxed_rank_001_*.pdb"))
        c08, c09 = contact_set(p08), contact_set(p09)
        union = c08 | c09
        pose.append({"complex_id": ident, "receptor_aligned_peptide_backbone_CA_RMSD_seed09_to_seed08_A": receptor_aligned_peptide_rmsd(p08, p09),
                     "contact_pair_jaccard_seed08_vs_seed09": len(c08 & c09) / len(union) if union else None,
                     "interpretation": "descriptive pose/geometry only; no threshold and no affinity inference"})

    execution_path = PILOT / "run_execution_record_v1.json"
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    failures = [e for e in execution["events"] if e.get("returncode") not in (0, None) or e.get("stop_reason")]
    execution.update({
        "status": "COMPLETED_WITH_STRUCTURAL_GATE_FAILURES_AND_ESM_LOADER_COMPAT_RETRY",
        "completed_at": execution.get("updated_at"),
        "gpu_wall_seconds": round(sum(float(e["duration_seconds"]) for e in execution["events"] if e.get("gpu_metered")), 3),
        "cpu_wall_seconds_recorded": round(sum(float(e["duration_seconds"]) for e in execution["events"] if not e.get("gpu_metered")), 3),
        "successful_prediction_seed_jobs": sum(e.get("label", "").startswith("afm_seed_") and e.get("returncode") == 0 for e in execution["events"]),
        "successful_complex_scoring_rows": {"GraphPep": len(graphpep), "FoldX": len(foldx), "ESM2_sequence_prior": len(esm)},
        "structural_quality_pass_cells": gate["quality_pass_cells"],
        "direction_evaluable_pairs": 0,
        "failed_events_preserved": failures,
        "final_output_bytes_delta": None,
    })
    pilot_bytes = sum(p.stat().st_size for p in PILOT.rglob("*") if p.is_file())
    execution["final_output_bytes_delta"] = pilot_bytes - int(execution["output_baseline_bytes"])
    execution["output_cap_pass"] = execution["final_output_bytes_delta"] <= execution["output_cap_bytes"]
    execution["gpu_budget_pass"] = execution["gpu_wall_seconds"] <= int(config["budgets"]["max_gpu_hours"] * 3600)
    execution_path.write_text(json.dumps(execution, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    summary = {
        "schema": "bclw_six_complex_pilot_paired_summary_v1",
        "status": "ENGINEERING_PILOT_COMPLETE_NO_VALID_COMPLEX_DIRECTION",
        "config_sha256": sha(PILOT / "run_config_frozen_v1.json"),
        "expected_cells": 6,
        "observed_cells": gate["observed_cells"],
        "identity_pass_cells": sum(bool(r.get("identity_pass")) for r in gate["rows"]),
        "contact_pass_cells": sum(bool(r.get("contact_pass")) for r in gate["rows"]),
        "quality_pass_cells": gate["quality_pass_cells"],
        "graphpep_success_rows": len(graphpep),
        "foldx_success_rows": sum(r["status"] == "ok" for r in foldx.values()),
        "esm_sequence_rows": len(esm),
        "gpu_wall_seconds_including_failed_initial_ESM_load": execution["gpu_wall_seconds"],
        "cpu_wall_seconds_for_relax_and_FoldX": execution["cpu_wall_seconds_recorded"],
        "new_output_bytes_delta": execution["final_output_bytes_delta"],
        "output_cap_bytes": execution["output_cap_bytes"],
        "gpu_budget_pass": execution["gpu_budget_pass"],
        "output_cap_pass": execution["output_cap_pass"],
        "failed_events_preserved": failures,
        "paired_proxy_diagnostics": paired,
        "sequence_prior_diagnostics": esm_pairs,
        "cross_seed_pose_diagnostics": pose,
        "interpretation_boundary": "All six complex cells failed the frozen AFM confidence gate, so none of the GraphPep/FoldX paired score differences is direction-evaluable under the frozen rule. Report raw scores and deltas only as diagnostics on failed structures. ESM2 PLL is a sequence-only prior. None of these outputs is experimental Kd or complete binding free energy.",
    }
    target = PILOT / "paired_proxy_summary_v1.json"
    target.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": summary["status"], "identity_pass": summary["identity_pass_cells"],
        "contact_pass": summary["contact_pass_cells"], "quality_pass": summary["quality_pass_cells"],
        "GraphPep_rows": summary["graphpep_success_rows"], "FoldX_rows": summary["foldx_success_rows"],
        "ESM_rows": summary["esm_sequence_rows"], "gpu_seconds": summary["gpu_wall_seconds_including_failed_initial_ESM_load"],
        "cpu_seconds": summary["cpu_wall_seconds_for_relax_and_FoldX"], "output_delta_bytes": summary["new_output_bytes_delta"],
        "paired_proxy_diagnostics": paired, "sequence_prior_diagnostics": esm_pairs, "cross_seed_pose": pose,
    }, ensure_ascii=False, indent=2))
    return 0 if summary["observed_cells"] == 6 and summary["output_cap_pass"] and summary["gpu_budget_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
