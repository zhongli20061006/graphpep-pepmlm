#!/usr/bin/env bash
# STEP 1 (full): generate 1000 candidate peptides for NCAM1 with PepMLM-650M.
set -u
PY=/home/g/miniconda3/envs/graphpep_env/bin/python
export TOKENIZERS_PARALLELISM=false
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
ROOT=/home/g/pepmlm_pipeline
mkdir -p $ROOT/logs $ROOT/output

echo "=================== STEP1 START $(date -Is) ==================="
echo "NOTE: output/ contains a corrupt 1.5GB truncated ESMFold download; removing it."
"$ROOT/scripts/sync_pipeline.sh" >/dev/null 2>&1 || true

cd $ROOT
$PY pipeline/step1_generate.py \
    --n 1000 --length 10 --top-k 3 --seed 42 \
    --dtype fp16 \
    --sample-batch 100 \
    --ppl-chunk 40 --ppl-row-chunk 10 \
    --verify-ppl 3 \
    --out $ROOT/output/generated_peptides.csv 2>&1 | grep -vE 'pkg_resources|declare_namespace'
echo "=================== STEP1 END $(date -Is) exit=$? ==================="
