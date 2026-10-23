#!/usr/bin/env bash
# Prove the patched resmap fixes GraphPep's -each path, using the same argv GraphPep.sh uses.
set -u
GP=/home/g/GraphPep_v1.3
W=/home/g/pepmlm_pipeline/work/gp_each_test
rm -rf "$W" && mkdir -p "$W"
cd "$W"

cp -f $GP/example/1B9J/protein.pdb pro_ori.pdb
cp -f $GP/example/1B9J/decoys.pdb  pep_ori.pdb
# renumbered/cleaned inputs are what GraphPep.sh feeds resmap; for this test the
# example files are already clean enough that renumming is identity-ish
cp -f pro_ori.pdb pro_renum.pdb
cp -f pep_ori.pdb pep_renum.pdb
cp -f $GP/example/1B9J/score_each.csv tmp_each.csv

echo "=== ORIGINAL resmap.py ==="
/home/g/miniconda3/envs/graphpep_env/bin/python $GP/bin/resmap.py \
    pro_ori.pdb pro_renum.pdb pep_ori.pdb pep_renum.pdb tmp_each.csv out_orig.csv 2>&1 | tail -4
echo "exit=$? ; out_orig.csv exists: $([ -f out_orig.csv ] && echo YES || echo NO)"
echo

echo "=== PATCHED resmap_compat.py ==="
/home/g/miniconda3/envs/graphpep_env/bin/python \
    /home/g/pepmlm_pipeline/pipeline/graphpep_patches/resmap_compat.py \
    pro_ori.pdb pro_renum.pdb pep_ori.pdb pep_renum.pdb tmp_each.csv out_patched.csv 2>&1 | tail -6
echo "out_patched.csv exists: $([ -f out_patched.csv ] && echo YES || echo NO)"
echo
if [ -f out_patched.csv ]; then
  echo "--- out_patched.csv (head) ---"
  head -6 out_patched.csv
  echo "rows: $(wc -l < out_patched.csv)"
fi
