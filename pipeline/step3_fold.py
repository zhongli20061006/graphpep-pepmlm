#!/usr/bin/env python3
"""Step 3a - fold each candidate peptide together with the receptor using ESMFold.

STATUS / SCOPE WARNING
----------------------
ESMFold is a single-chain structure predictor. It was not trained on protein-peptide
interfaces and outputs no interface-confidence signal, so the relative placement of the
peptide against the receptor is not a reliable binding pose. This script therefore does
NOT claim to produce bound complexes. It is provided so that, once the ESMFold weights are
available, the protocol can be run AND audited:

  * every complex is folded with the receptor at the N-terminus and the peptide at the
    C-terminus of a single chain, so both parts come from one consistent coordinate frame;
  * for each complex we measure the minimum inter-chain distance and the number of
    peptide residues within GraphPep's 6 A contact threshold (bin/config.py);
  * complexes with no contact produce no interaction graph, and their GraphPep score is
    meaningless - this script reports the zero-contact rate so that is never hidden.

The receptor is folded once and cached, since it is identical for every peptide.

Outputs:
  output/pdbs/{peptide_id}_{receptor,peptide,complex}.pdb
  output/fold_contacts.csv      (peptide_id, min_dist, n_contact_res, contact_ok)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch

ROOT = "/home/g/pepmlm_pipeline"
ESMFOLD_DIR = os.path.join(ROOT, "models", "esmfold_v1")
PDB_DIR = os.path.join(ROOT, "output", "pdbs")
CACHE_DIR = os.path.join(ROOT, "cache", "fold")
ERR_LOG = os.path.join(ROOT, "output", "error.log")

CONTACT_CUTOFF = 6.0      # GraphPep's dis_threshold (bin/config.py)


def log_err(msg: str) -> None:
    os.makedirs(os.path.dirname(ERR_LOG), exist_ok=True)
    with open(ERR_LOG, "a") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
    print(f"ERROR: {msg}", file=sys.stderr)


def read_fasta(path: str) -> str:
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def write_pdb(out_path: str, coords: np.ndarray, seq: str, chain: str, start_res: int = 1):
    """Minimal backbone-only PDB writer (N, CA, C, O)."""
    three = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
        "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
        "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
        "Y": "TYR", "V": "VAL",
    }
    atoms = ["N", "CA", "C", "O"]
    lines = []
    serial = 1
    for i, aa in enumerate(seq):
        resn = three.get(aa, "UNK")
        for ai, atom in enumerate(atoms):
            x, y, z = coords[i, ai]
            lines.append(
                f"ATOM  {serial:5d} {atom:<4s}{resn:>3s} {chain}{start_res+i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           {atom[0]:>2s}")
            serial += 1
    lines.append("END")
    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def extract_backbone(pdb_text: str) -> tuple[np.ndarray, str]:
    """Return (N, L, 4, 3) backbone coords and the one-letter sequence from an ESMFold PDB."""
    order = ["N", "CA", "C", "O"]
    res_coords: dict[int, dict[str, tuple[float, float, float]]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom = line[12:16].strip()
        if atom not in order:
            continue
        resseq = int(line[22:26])
        res_coords.setdefault(resseq, {})[atom] = (
            float(line[30:38]), float(line[38:46]), float(line[46:54]))
    keys = sorted(res_coords)
    coords = np.zeros((len(keys), 4, 3), dtype=np.float64)
    for i, k in enumerate(keys):
        for ai, atom in enumerate(order):
            coords[i, ai] = res_coords[k][atom]
    return coords, "".join("X" for _ in keys)


def min_inter_chain_distance(rec: np.ndarray, pep: np.ndarray) -> tuple[float, int]:
    """Minimum all-atom distance, and the number of peptide residues within CONTACT_CUTOFF."""
    r = rec.reshape(-1, 3)
    p = pep.reshape(-1, 3)
    d = np.linalg.norm(r[:, None, :] - p[None, :, :], axis=-1)
    min_d = float(d.min())
    per_res = d.min(axis=0).reshape(pep.shape[0], -1).min(axis=1)
    n_contact = int((per_res <= CONTACT_CUTOFF).sum())
    return min_d, n_contact


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default=os.path.join(ROOT, "output", "top20_lists_combined.csv"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--receptor-range", type=str, default="",
                    help="1-based inclusive 'start-end' of the receptor construct, e.g. 20-250")
    args = ap.parse_args()

    if not os.path.isfile(os.path.join(ESMFOLD_DIR, "pytorch_model.bin")):
        msg = (f"ESMFold weights missing at {ESMFOLD_DIR}/pytorch_model.bin - stage 3 cannot "
               f"run. Download facebook/esmfold_v1 (8.44 GB) first.")
        log_err(msg)
        print(msg)
        return 2

    os.makedirs(PDB_DIR, exist_ok=True)
    os.makedirs(CACHE_DIR, exist_ok=True)
    t0 = time.time()

    print("loading ESMFold ...", flush=True)
    import esm

    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()
    model.set_chunk_size(64)          # keeps the folding trunk inside 8 GB

    target = read_fasta(os.path.join(ROOT, "data/targets/NCAM1.fasta"))
    if args.receptor_range:
        a, b = (int(x) for x in args.receptor_range.split("-"))
        receptor = target[a - 1:b]
    else:
        receptor = target
    print(f"  receptor construct: {len(receptor)} aa (target full length {len(target)})")

    pep_df = pd.read_csv(args.scores) if os.path.isfile(args.scores) else None
    if pep_df is None:
        msg = f"missing peptide list: {args.scores}"
        log_err(msg)
        print(msg)
        return 2
    if args.limit:
        pep_df = pep_df.head(args.limit)

    # peptide-only folds are cheap and cached; the receptor is folded once
    rows = []
    for _, r in pep_df.iterrows():
        pid, pep = r["peptide_id"], r["peptide_sequence"]
        chain = receptor + pep
        cache = os.path.join(CACHE_DIR, f"{pid}.pdb")
        if not os.path.isfile(cache):
            with torch.no_grad():
                out = model.infer_pdb(chain)
            with open(cache, "w") as fh:
                fh.write(out)
            torch.cuda.empty_cache()
        coords, _ = extract_backbone(open(cache).read())
        n_rec = len(receptor)
        rec_c = coords[:n_rec]
        pep_c = coords[n_rec:n_rec + len(pep)]
        if len(pep_c) != len(pep):
            log_err(f"{pid}: folded peptide length {len(pep_c)} != {len(pep)}; skipping")
            continue
        write_pdb(os.path.join(PDB_DIR, f"{pid}_receptor.pdb"), rec_c, receptor, "A")
        write_pdb(os.path.join(PDB_DIR, f"{pid}_peptide.pdb"), pep_c, pep, "B")
        # combined complex, two chains, so GraphPep sees a real two-chain input
        with open(os.path.join(PDB_DIR, f"{pid}_complex.pdb"), "w") as fh:
            fh.write(open(os.path.join(PDB_DIR, f"{pid}_receptor.pdb")).read().replace("END\n", ""))
            fh.write(open(os.path.join(PDB_DIR, f"{pid}_peptide.pdb")).read())
        min_d, n_contact = min_inter_chain_distance(rec_c, pep_c)
        rows.append({"peptide_id": pid, "source": r.get("source", ""),
                     "peptide_sequence": pep, "min_dist": min_d,
                     "n_contact_res": n_contact,
                     "contact_ok": n_contact > 0})
        print(f"  {pid} {pep}  min_dist={min_d:6.2f} A  contact_res={n_contact}/{len(pep)}")

    df = pd.DataFrame(rows)
    out = os.path.join(ROOT, "output", "fold_contacts.csv")
    df.to_csv(out, index=False)
    n_zero = int((~df["contact_ok"]).sum()) if len(df) else 0
    print(f"\nSTEP3a DONE time={time.time()-t0:.1f}s")
    print(f"  folded {len(df)} complexes -> {PDB_DIR}")
    print(f"  ZERO-CONTACT complexes: {n_zero}/{len(df)} "
          f"({100*n_zero/max(len(df),1):.1f}%)  <-- these have uninformative GraphPep scores")
    print(f"  contact table -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
