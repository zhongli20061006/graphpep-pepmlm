#!/usr/bin/env bash
# Verify step 1 output before trusting it.
set -u
F=/home/g/pepmlm_pipeline/output/generated_peptides.csv
echo "=== file ==="
ls -la "$F"
echo
echo "=== header + first rows ==="
head -5 "$F"
echo "..."
echo
/home/g/miniconda3/envs/graphpep_env/bin/python - <<'PY'
import pandas as pd
f = "/home/g/pepmlm_pipeline/output/generated_peptides.csv"
df = pd.read_csv(f)
print("=== integrity checks ===")
print("  rows                 :", len(df))
print("  columns              :", list(df.columns))
print("  unique peptide_id    :", df.peptide_id.nunique())
print("  unique sequences     :", df.peptide_sequence.nunique())
print("  duplicate sequences  :", int(df.peptide_sequence.duplicated().sum()))
lens = df.peptide_sequence.str.len().unique().tolist()
print("  sequence lengths     :", lens)
print("  null in any column   :", int(df.isna().sum().sum()))
# only the 20 canonical amino acids?
import re
bad = sorted({c for s in df.peptide_sequence for c in s if c not in "ACDEFGHIKLMNPQRSTVWY"})
print("  non-canonical chars  :", bad if bad else "none")
print()
print("  ppl_score stats:")
print(df.ppl_score.describe().to_string())
print()
print("=== amino-acid composition (top 8) ===")
from collections import Counter
c = Counter("".join(df.peptide_sequence))
tot = sum(c.values())
for aa, n in c.most_common(8):
    print(f"  {aa}: {100*n/tot:5.1f}%")
print()
print("=== top 10 by PPL (lowest) ===")
print(df.nsmallest(10, "ppl_score").to_string(index=False))
PY
