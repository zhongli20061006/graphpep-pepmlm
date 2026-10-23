"""Read-only independent inspection of previously completed structures and terms."""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.vectors import Vector, calc_dihedral

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
PARSER = PDBParser(QUIET=True)


def load(pid, seed, model):
    path = ROOT / "output/foldx_affinity_mdm2_20261006/jobs" / f"{pid}_seed{seed}_model{model}" / "input_Repair.pdb"
    structure = PARSER.get_structure(path.stem, str(path))[0]
    assert [len(list(structure[c].get_residues())) for c in ("A", "B")] == [85, 12]
    return structure


def torsions(chain):
    result = []
    for n in range(1, 13):
        r = chain[n]
        phi = math.degrees(calc_dihedral(Vector(chain[n - 1]["C"].coord), Vector(r["N"].coord), Vector(r["CA"].coord), Vector(r["C"].coord))) if n > 1 else None
        psi = math.degrees(calc_dihedral(Vector(r["N"].coord), Vector(r["CA"].coord), Vector(r["C"].coord), Vector(chain[n + 1]["N"].coord))) if n < 12 else None
        result.append(dict(position=n, residue=r.resname, phi_deg=phi, psi_deg=psi, alpha_like=phi is not None and psi is not None and -100 <= phi <= -30 and -80 <= psi <= -5))
    return result


def aligned_peptide_rmsd(wt, mutant):
    fixed = np.array([wt["A"][n]["CA"].coord for n in range(1, 86)], dtype=float)
    moving = np.array([mutant["A"][n]["CA"].coord for n in range(1, 86)], dtype=float)
    fc, mc = fixed.mean(0), moving.mean(0)
    u, _, vh = np.linalg.svd((moving - mc).T @ (fixed - fc))
    if np.linalg.det(u @ vh) < 0:
        u[:, -1] *= -1
    rotation = u @ vh
    difference = np.array([(mutant["B"][n]["CA"].coord - mc) @ rotation + fc - wt["B"][n]["CA"].coord for n in range(1, 13)])
    result = {"all_CA_rmsd_A": float(np.sqrt(np.mean(np.sum(difference ** 2, axis=1)))), "B8_12_CA_rmsd_A": float(np.sqrt(np.mean(np.sum(difference[7:] ** 2, axis=1))))}
    for prefix, chain, first in (("A_backbone", "A", 1), ("B_backbone", "B", 1), ("B8_12_backbone", "B", 8)):
        keys = [(r.id[1], a.name) for r in wt[chain] if r.id[1] >= first for a in r if a.name in ("N", "CA", "C", "O", "OXT")]
        delta = np.array([(mutant[chain][n][a].coord - mc) @ rotation + fc - wt[chain][n][a].coord for n, a in keys])
        result[prefix + "_rmsd_A"] = float(np.sqrt(np.mean(np.sum(delta ** 2, axis=1))))
    return result


def main():
    protocol = json.loads((OUT / "protocol.json").read_text())
    for name, value in protocol["source_hashes"].items():
        with (ROOT / name).open("rb") as stream:
            assert hashlib.file_digest(stream, "sha256").hexdigest() == value
    structures, comparisons = [], []
    for seed in (20261011, 20261012):
        for model in (1, 2, 3):
            wt = load("P53_WT", seed, model)
            for pid in ("P53_WT", "P53_P27A", "P53_L22A"):
                s = wt if pid == "P53_WT" else load(pid, seed, model)
                angles = torsions(s["B"])
                structures.append(dict(peptide_id=pid, seed=seed, model=model, torsions=angles,
                                       alpha_like_B8_11=sum(r["alpha_like"] for r in angles[7:11]),
                                       O7_N11_distance_A=float(np.linalg.norm(s["B"][7]["O"].coord - s["B"][11]["N"].coord)),
                                       N11_is_possible_amide_donor=s["B"][11].resname != "PRO"))
                if pid != "P53_WT":
                    comparisons.append(dict(peptide_id=pid, seed=seed, model=model, **aligned_peptide_rmsd(wt, s)))
    with (ROOT / protocol["source_energy_terms_csv"]).open(newline="") as stream:
        terms = {(r["job_id"], r["repeat"], r["role"], r["term"]): float(r["value"]) for r in csv.DictReader(stream) if Path(r["table_path"]).name.startswith("Interaction_")}
    with (ROOT / protocol["source_pair_csv"]).open(newline="") as stream:
        pairs = [r for r in csv.DictReader(stream) if r["peptide_id"] in ("P53_P27A", "P53_L22A")]
    differences = {}
    for r in pairs:
        assert r["status"] == "ok"
        keys = [key[3] for key in terms if key[:3] == (r["job_id"], r["repeat"], "wt_reference")]
        for term in keys:
            value = terms[r["job_id"], r["repeat"], "wt_reference", term] - terms[r["job_id"], r["repeat"], "mutant", term]
            differences.setdefault(r["peptide_id"], {}).setdefault(term, []).append(value)
    summaries = {pid: {term: statistics.median(values) for term, values in records.items()} for pid, records in differences.items()}
    result = dict(structures=structures, aligned_comparisons=comparisons, pooled_descriptive_term_contrasts=summaries,
                  caution="Geometry proxies and pooled component medians are exploratory; not measured helicity, no energy-component sum assumption, not the primary aggregation")
    (OUT / "existing_structure_inspection.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"structures": len(structures), "comparisons": comparisons,
                      "P27A_largest_component_medians": sorted(summaries["P53_P27A"].items(), key=lambda item: -abs(item[1]))[:7]}, indent=2))


if __name__ == "__main__":
    main()
