#!/usr/bin/env bash
set -u
PY=/home/g/miniconda3/envs/graphpep_env/bin/python
export TOKENIZERS_PARALLELISM=false
cd /home/g/pepmlm_pipeline
rm -f cache/ppl_cache.json
echo "=== fp16 validation run (40 peptides) $(date +%H:%M:%S) ==="
time $PY pipeline/step1_generate.py \
    --n 40 --length 10 --top-k 3 --dtype fp16 \
    --sample-batch 40 --ppl-chunk 32 --ppl-row-chunk 8 \
    --verify-ppl 3 \
    --out /home/g/pepmlm_pipeline/output/smoke_fp16.csv 2>&1 | \
  grep -vE 'pkg_resources|declare_namespace' | tail -25
