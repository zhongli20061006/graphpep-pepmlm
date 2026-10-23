#!/usr/bin/env bash
export HF_ENDPOINT=https://hf-mirror.com
OUT=/home/g/pepmlm_pipeline/probe
mkdir -p $OUT
for f in bindevaluator.py moppit.py scripts/bindevaluator.sh README.md; do
  echo "########## $f ##########"
  timeout 60 curl -4 -sSL "https://hf-mirror.com/ChatterjeeLab/moPPIt/resolve/main/$f" -o "$OUT/$(basename $f)" 2>&1
  wc -l "$OUT/$(basename $f)" 2>&1
done
echo "===== bindevaluator.py (head 150) ====="
head -150 $OUT/bindevaluator.py 2>&1
echo
echo "===== scripts/bindevaluator.sh ====="
cat $OUT/bindevaluator.sh 2>&1 | head -60
echo
echo "===== README.md ====="
head -80 $OUT/README.md 2>&1
