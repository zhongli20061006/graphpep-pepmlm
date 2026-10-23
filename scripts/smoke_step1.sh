#!/usr/bin/env bash
# Smoke test: tiny PepMLM run to validate the generator end-to-end and measure throughput.
set -u
PY=/home/g/miniconda3/envs/graphpep_env/bin/python
export TOKENIZERS_PARALLELISM=false
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True

echo "=== smoke test $(date -Is) ==="
nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader
echo
cd /home/g/pepmlm_pipeline
/usr/bin/time -v $PY pipeline/step1_generate.py \
    --n 20 --length 10 --top-k 3 \
    --sample-batch 20 --ppl-batch 10 \
    --out /home/g/pepmlm_pipeline/output/smoke_peptides.csv 2>&1 | \
  grep -vE '^\s+(Command being timed|User time|System time|Percent|Average|Major|Minor|Voluntary|Involuntary|Swaps|File system|Socket|Signals|Page size|Exit status|Elapsed \(wall)' | \
  head -60
echo
echo "=== peak GPU memory ==="
nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader
echo
echo "=== output ==="
cat /home/g/pepmlm_pipeline/output/smoke_peptides.csv 2>&1 | head -10
