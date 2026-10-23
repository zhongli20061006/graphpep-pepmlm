#!/usr/bin/env python3
"""Reproduce the resmap.py KeyError(0) precisely, and test the minimal fix.

resmap.py does:
    df = pd.read_csv(csv_in, dtype=str)      # header is 'DECOY_ID,RESID_PEP,RESID_PRO,CONFIDENT_SCORES'
    df.iloc[:,1] = df.iloc[:,1].str.replace(...)
    for _, row in df.iterrows():
        decoy = row[0]        # <-- KeyError: 0

The question is why row[0] fails here but works on a toy frame in pandas 3.0.5.
"""
import pandas as pd

print("pandas", pd.__version__)

# Load the real tmp_each.csv that GraphPep writes (the model's per-residue output)
import glob
cands = glob.glob("/home/g/pepmlm_pipeline/work/gp_test/tmp_each*") + \
        glob.glob("/tmp/tmp.*/tmp_each.csv")
print("candidate tmp_each files:", cands)

path = "/home/g/GraphPep_v1.3/example/1B9J/score_each.csv"
print(f"\n--- trying the repo's stored per-residue file: {path} ---")
df = pd.read_csv(path, dtype=str)
print("columns:", list(df.columns))
print("index:", df.index[:4].tolist(), "dtype:", df.index.dtype)
row = next(df.iterrows())[1]
print("row.index:", list(row.index))
try:
    print("row[0] ->", row[0])
except Exception as e:
    print(f"row[0] FAILS -> {type(e).__name__}: {e}")
try:
    print("row.iloc[0] ->", row.iloc[0])
except Exception as e:
    print(f"row.iloc[0] FAILS -> {type(e).__name__}: {e}")

# Now WITHOUT a header (which is what the injected tmp_each.csv looks like when the
# '_outputrescsv' writer appends without a header on non-first batches): the first data
# row is consumed as the header, so columns become ['1','1:B','32:A','0.9993']
print("\n--- simulate headerless/odd-header read ---")
df2 = pd.read_csv(path, dtype=str, header=None)
print("columns:", list(df2.columns), "shape:", df2.shape)
row2 = next(df2.iterrows())[1]
try:
    print("row2[0] ->", row2[0])
except Exception as e:
    print(f"row2[0] FAILS -> {type(e).__name__}: {e}")
print("row2.iloc[0] ->", row2.iloc[0])
