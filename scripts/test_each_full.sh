#!/usr/bin/env bash
# Prove the corrected resmap handles the real "resid:chain" format, then run the FULL
# GraphPep -each pipeline with the patched resmap substituted in.
set -u
GP=/home/g/GraphPep_v1.3
PATCH=/home/g/pepmlm_pipeline/pipeline/graphpep_patches/resmap_compat.py
PY=/home/g/miniconda3/envs/graphpep_env/bin/python

W=/home/g/pepmlm_pipeline/work/gp_each_full
rm -rf "$W" && mkdir -p "$W"
cd "$W"
cp -f $GP/example/1B9J/protein.pdb .
cp -f $GP/example/1B9J/decoys.pdb .

echo "=== 1. corrected resmap on the authors' own per-residue output ==="
cp -f $GP/example/1B9J/score_each.csv tmp_each.csv
$PY $PATCH protein.pdb protein.pdb decoys.pdb decoys.pdb tmp_each.csv out_each.csv 2>&1 | tail -6
echo "--- out_each.csv (head) ---"
head -6 out_each.csv 2>&1
echo

echo "=== 2. FULL GraphPep -each via a patched copy of GraphPep.sh ==="
# Make a pipeline-owned copy of GraphPep.sh that calls our resmap instead of the broken one
sed "s#\$GraphPep_root/bin/resmap.py#$PATCH#g" $GP/GraphPep.sh > gp_patched.sh
grep -n "resmap" gp_patched.sh
echo
timeout 1200 bash gp_patched.sh protein.pdb decoys.pdb \
    -out score_full.csv -each -outeach score_each_full.csv 2>&1 | tail -8
echo
echo "=== 3. results ==="
ls -la "$W"
echo "--- score_full.csv ---"; cat score_full.csv 2>&1
echo "--- score_each_full.csv (head) ---"; head -6 score_each_full.csv 2>&1
echo "--- rows: $(wc -l < score_each_full.csv 2>/dev/null || echo 0) ---"
echo
echo "=== 4. compare simple score to reference ==="
echo "reference: $(tr '\n' ' ' < $GP/example/1B9J/score.csv)"
echo "ours:      $(tr '\n' ' ' < score_full.csv 2>/dev/null)"
