#!/usr/bin/env python3
"""Characterise the generated candidate pool (diversity + composition bias).

Saved to output/pool_stats.json and folded into the report so the reader can judge how
diverse the PepMLM output actually is before trusting any downstream comparison.
"""
from __future__ import annotations

import json
import os
from collections import Counter

import numpy as np
import pandas as pd

ROOT = "/home/g/pepmlm_pipeline"
OUT = os.path.join(ROOT, "output")
AA = "ACDEFGHIKLMNPQRSTVWY"


def main() -> int:
    df = pd.read_csv(os.path.join(OUT, "generated_peptides.csv"))
    seqs = df["peptide_sequence"].tolist()
    n, nuniq = len(df), df["peptide_sequence"].nunique()

    stats: dict = {"n_candidates": n, "n_unique": nuniq,
                   "duplicate_fraction": float(1 - nuniq / n)}

    # per-position amino-acid entropy, normalised by log(20)
    L = len(seqs[0])
    ent, cons = [], []
    for pos in range(L):
        c = Counter(s[pos] for s in seqs)
        p = np.array([c.get(aa, 0) for aa in AA], dtype=float)
        p = p / p.sum()
        nz = p[p > 0]
        ent.append(float(-(nz * np.log(nz)).sum() / np.log(len(AA))))
        cons.append(c.most_common(1)[0])
    stats["mean_positional_entropy"] = float(np.mean(ent))
    stats["positional_entropy"] = [round(e, 4) for e in ent]
    stats["position_consensus"] = [f"{aa}({100*k/n:.0f}%)" for aa, k in cons]

    # overall composition
    comp = Counter("".join(seqs))
    tot = sum(comp.values())
    stats["composition_pct"] = {aa: round(100 * comp.get(aa, 0) / tot, 2) for aa in AA}

    # diversity of the top-N by PPL vs the whole pool
    for tag, sub in (("all", df), ("top20_ppl", df.nsmallest(20, "ppl_score")),
                     ("top100_ppl", df.nsmallest(100, "ppl_score"))):
        s = sub["peptide_sequence"].tolist()
        cc = Counter("".join(s))
        tt = sum(cc.values())
        stats[f"composition_pct_{tag}"] = {aa: round(100 * cc.get(aa, 0) / tt, 2) for aa in AA}
        stats[f"unique_{tag}"] = int(sub["peptide_sequence"].nunique())

    with open(os.path.join(OUT, "pool_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=2)

    print("=== candidate pool ===")
    print(f"  candidates={n}  unique={nuniq}  duplicate fraction={stats['duplicate_fraction']:.3f}")
    print(f"  mean positional entropy (0-1, 1=diverse): {stats['mean_positional_entropy']:.3f}")
    print(f"  per-position entropy: {stats['positional_entropy']}")
    print(f"  per-position consensus: {stats['position_consensus']}")
    print()
    print("=== overall composition (top 10) ===")
    top = sorted(stats["composition_pct"].items(), key=lambda kv: -kv[1])[:10]
    for aa, pct in top:
        print(f"  {aa}: {pct:5.2f}%")
    print()
    print("=== composition of top100 by PPL (top 10) ===")
    top2 = sorted(stats["composition_pct_top100_ppl"].items(), key=lambda kv: -kv[1])[:10]
    for aa, pct in top2:
        print(f"  {aa}: {pct:5.2f}%")
    print()
    print(f"written: {os.path.join(OUT, 'pool_stats.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
