#!/usr/bin/env bash
# Quantify the GraphPep pandas-3 incompatibility:
#   1. does the simple (non -each) path still work?  -> gives PRED_SCORE only
#   2. exactly which lines in resmap.py break, and what pandas versions allow row[0]
set -u
GP=/home/g/GraphPep_v1.3
CDIR=$GP/example/1B9J
WORK=/tmp/gp_validate2
rm -rf $WORK && mkdir -p $WORK
cd "$CDIR"

echo "=== A. simple path (no -each): PRED_SCORE only ==="
timeout 900 bash $GP/GraphPep.sh protein.pdb decoys.pdb -out $WORK/score_simple.csv 2>&1 | tail -6
echo "--- result ---"
cat $WORK/score_simple.csv 2>&1
echo
echo "repo reference:"; cat $CDIR/score.csv
echo
echo "=== B. resmap.py positional-index usages (pandas 3 removed row[0] on Series) ==="
grep -n "row\[" $GP/bin/resmap.py
echo
echo "Pandas version in graphpep_env:"
/home/g/miniconda3/envs/graphpep_env/bin/python -c "import pandas; print('  pandas', pandas.__version__)"
echo
echo "=== C. confirm the failure mode in isolation ==="
/home/g/miniconda3/envs/graphpep_env/bin/python - <<'PY'
import pandas as pd
print("pandas", pd.__version__)
df = pd.DataFrame([["1","2","3","0.5"]])
row = df.iloc[0]
try:
    print("row[0] ->", row[0])
except Exception as e:
    print(f"row[0] FAILS: {type(e).__name__}: {e}")
try:
    print("row.iloc[0] ->", row.iloc[0])
except Exception as e:
    print("row.iloc[0] fails too:", e)
PY
echo
echo "=== D. GraphPep's declared pandas pin ==="
grep -riE "pandas" $GP/requirements.txt $GP/environment.yml 2>/dev/null
echo
echo "=== E. is a pandas<2 install available in the mirror (dry-run resolve) ==="
timeout 300 /home/g/miniconda3/envs/graphpep_env/bin/pip index versions pandas 2>&1 | head -5
