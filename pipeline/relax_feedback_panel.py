"""Apply identical restrained Amber minimization to every AFM complex."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

os.environ.setdefault("JAX_PLATFORMS", "cpu")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("OPENMM_CPU_THREADS", "2")

import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from colabfold.relax import relax_me

PARSER = PDBParser(QUIET=True)
SETTINGS = {"max_iterations": 200, "tolerance": 2.39, "stiffness": 10.,
            "max_outer_iterations": 3, "use_gpu": False}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry(path: Path) -> tuple[dict, np.ndarray]:
    model = PARSER.get_structure(path.stem, str(path))[0]
    if set(chain.id for chain in model) != {"A", "B"}:
        raise ValueError(f"Unexpected chain IDs in {path}")
    sequences, heavy, cas = {}, {}, []
    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]
        sequences[chain.id] = "".join(seq1(r.resname) for r in residues)
        cas.extend(r["CA"].coord for r in residues)
        heavy[chain.id] = np.asarray([atom.coord for r in residues for atom in r
                                     if atom.element.upper() not in {"H", "D"}], dtype=float)
    distances = np.linalg.norm(heavy["A"][:, None] - heavy["B"][None, :], axis=2)
    return {"sequences": sequences, "minimum_heavyatom_distance": float(distances.min()),
            "clashing_atom_pairs_lt_1p5A": int((distances < 1.5).sum())}, np.asarray(cas)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    files = sorted(args.input.glob("*_unrelaxed_rank_*.pdb"))
    if not files:
        raise ValueError("No AFM unrelaxed structures found")
    manifest_path = args.output / "relaxation_manifest.json"
    prior = json.loads(manifest_path.read_text())["models"] if manifest_path.exists() else []
    cache = {row["source_path"]: row for row in prior}
    rows = []
    for index, source in enumerate(files, 1):
        destination = args.output / source.name.replace("_unrelaxed_", "_relaxed_")
        source_hash = sha256(source)
        saved = cache.get(str(source), {})
        if (destination.exists() and saved.get("source_sha256") == source_hash and
                saved.get("relaxed_sha256") == sha256(destination) and saved.get("settings") == SETTINGS):
            rows.append(saved)
            print(f"Amber {index}/{len(files)} cached: {source.name}", flush=True)
            continue
        started = time.perf_counter()
        before, before_ca = geometry(source)
        destination.write_text(relax_me(pdb_filename=str(source), **SETTINGS))
        after, after_ca = geometry(destination)
        if before["sequences"] != after["sequences"] or before_ca.shape != after_ca.shape:
            raise ValueError("Amber changed chain identities or residue count")
        confidence = source.with_name(source.name.replace("_unrelaxed_", "_scores_")).with_suffix(".json")
        shutil.copy2(confidence, args.output / confidence.name)
        row = {"source_path": str(source), "relaxed_path": str(destination),
               "source_sha256": source_hash, "relaxed_sha256": sha256(destination),
               "confidence_sha256_unchanged": sha256(confidence), "settings": SETTINGS,
               "before": before, "after": after,
               "all_chain_CA_displacement_RMS_A": float(np.sqrt(np.mean(np.sum((after_ca - before_ca) ** 2, axis=1)))),
               "seconds": time.perf_counter() - started}
        rows.append(row)
        manifest_path.write_text(json.dumps({"settings": SETTINGS, "models": rows}, indent=2) + "\n")
        print(f"Amber {index}/{len(files)}: min-distance {before['minimum_heavyatom_distance']:.3f} -> {after['minimum_heavyatom_distance']:.3f} A, CA displacement {row['all_chain_CA_displacement_RMS_A']:.3f} A", flush=True)
    print(f"Relaxed {len(rows)} complexes", flush=True)


if __name__ == "__main__":
    main()
