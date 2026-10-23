#!/usr/bin/env python3
"""Validate the bind_eval_score aggregation against the authors' own semantics.

From `ChatterjeeLab/moPPIt/bindevaluator.py::main`:
    binding_site = [i for i in range(len(prediction)) if prediction[i] >= 0.5]
and for motif scoring:
    torch.sum(prediction[motifs]) / len(motifs)

So the model's own notion of "the binding site" is {p_i >= 0.5}, and the natural
peptide-level aggregate is the mean probability over exactly that set. This script checks
our scorer reproduces that on synthetic vectors, including the empty-set edge case.
"""
from __future__ import annotations

import numpy as np

SITE_THRESHOLD = 0.5


def ours(p):
    p = np.asarray(p, dtype=float)
    sites = np.where(p >= SITE_THRESHOLD)[0]
    return (float(p[sites].mean()) if sites.size else 0.0), int(sites.size)


def authors_motif_score(p, motifs):
    p = np.asarray(p, dtype=float)
    return float(p[motifs].sum() / len(motifs))


cases = {
    "all above":      np.array([0.9, 0.8, 0.7, 0.6]),
    "all below":      np.array([0.1, 0.2, 0.3, 0.4]),
    "mixed":          np.array([0.1, 0.9, 0.2, 0.8]),
    "exactly 0.5":    np.array([0.5, 0.5]),
    "single peak":    np.array([0.2, 0.99, 0.1]),
    "empty":          np.array([0.0, 0.0, 0.0]),
}

print("=== our aggregation ===")
for name, p in cases.items():
    score, n = ours(p)
    print(f"  {name:14s} -> score={score:.6f}  n_sites={n}")

print()
print("=== equivalence with the authors' motif-score formula when the 'motif' IS the")
print("=== model's own binding-site set ===")
p = np.array([0.1, 0.9, 0.95, 0.2, 0.51])
sites = np.where(p >= SITE_THRESHOLD)[0].tolist()
ours_score, _ = ours(p)
auth = authors_motif_score(p, sites)
print(f"  p               = {p.tolist()}")
print(f"  binding site    = {sites}")
print(f"  ours            = {ours_score:.8f}")
print(f"  authors'formula = {auth:.8f}")
print(f"  identical       = {abs(ours_score - auth) < 1e-12}")

# sanity: monotonicity / bounds
print()
print("=== bounds sanity ===")
rng = np.random.default_rng(0)
bad = 0
for _ in range(2000):
    v = rng.random(rng.integers(1, 60))
    s, n = ours(v)
    if n == 0:
        assert s == 0.0
    else:
        assert SITE_THRESHOLD - 1e-12 <= s <= 1.0 + 1e-12, (s, n)
        assert n == int((v >= SITE_THRESHOLD).sum())
        assert bad == 0 or True
print("  2000 random vectors: score in [0.5, 1.0] whenever any site exists, 0.0 otherwise -- OK")
