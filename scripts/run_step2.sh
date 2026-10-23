#!/usr/bin/env bash
# Run STEP 2 (BindEvaluator) once the checkpoint is complete.
# Usage: run_step2.sh [--limit N] [--peptides FILE]
set -u
PY=/home/g/miniconda3/envs/graphpep_env/bin/python
export TOKENIZERS_PARALLELISM=false
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
ROOT=/home/g/pepmlm_pipeline
CKPT=$ROOT/models/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt
EXPECT=2832618887

bash $ROOT/scripts/sync_pipeline.sh >/dev/null 2>&1 || true

if [ ! -f "$CKPT" ]; then
  echo "BindEvaluator checkpoint NOT READY at $CKPT"
  ls -la $CKPT.part 2>/dev/null
  exit 2
fi
SZ=$(stat -c %s "$CKPT")
echo "checkpoint size = $SZ (expected $EXPECT)"
if [ "$SZ" -ne "$EXPECT" ]; then
  echo "SIZE MISMATCH -> refusing to run"
  exit 3
fi

PEP=$ROOT/output/generated_peptides.csv
LIMIT=""
if [ "${1:-}" = "--limit" ]; then LIMIT="--limit $2"; fi

echo "=== STEP2 START $(date -Is) $LIMIT ==="
cd $ROOT
$PY pipeline/step2_bindeval.py \
    --peptides "$PEP" --out "$PEP" \
    $LIMIT \
    2>&1 | grep -vE 'pkg_resources|declare_namespace' | grep -vE 'UserWarning|warnings.warn'
echo "=== STEP2 END $(date -Is) exit=$? ==="
