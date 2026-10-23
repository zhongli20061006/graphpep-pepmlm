#!/usr/bin/env bash
# Smoke test 2: validate the batched pseudo-PPL path against the authors' reference loop,
# and measure the throughput improvement.
set -u
PY=/home/g/miniconda3/envs/graphpep_env/bin/python
export TOKENIZERS_PARALLELISM=false

echo "=== smoke2 $(date -Is) ==="
cd /home/g/pepmlm_pipeline
time $PY pipeline/step1_generate.py \
    --n 40 --length 10 --top-k 3 \
    --sample-batch 40 --ppl-chunk 32 --ppl-row-chunk 8 \
    --verify-ppl 3 \
    --out /home/g/pepmlm_pipeline/output/smoke2_peptides.csv 2>&1 | \
  grep -vE 'pkg_resources|declare_namespace' | tail -40
echo
echo "=== rows ==="
wc -l /home/g/pepmlm_pipeline/output/smoke2_peptides.csv
