#!/usr/bin/env bash
# Transport the Windows-side pipeline/*.py into WSL with LF endings, then syntax-check.
set -u
SRC=/mnt/e/Users/graphpep与pepmlm/pipeline
DST=/home/g/pepmlm_pipeline/pipeline
mkdir -p "$DST"
shopt -s nullglob
for f in "$SRC"/*.py; do
  base=$(basename "$f")
  tr -d '\r' < "$f" > "$DST/$base"
  echo "  -> $base ($(stat -c %s "$DST/$base") bytes)"
done
echo
echo "=== syntax check ==="
PY=/home/g/miniconda3/envs/graphpep_env/bin/python
for f in "$DST"/*.py; do
  if $PY -m py_compile "$f" 2>/tmp/pyerr.txt; then
    echo "  OK   $(basename "$f")"
  else
    echo "  FAIL $(basename "$f")"; cat /tmp/pyerr.txt
  fi
done
echo
echo "=== checks ==="
some_fail=0
check() {
  # 1 = must be present and non-empty, 0 = may be absent
  local path="$1" must="$2" label="$3"
  if [ -s "$path" ]; then
    echo "  PRESENT  $label  ($(stat -c %s "$path") bytes)"
  elif [ "$must" = "1" ]; then
    echo "  MISSING  $label   <-- required"; some_fail=1
  else
    echo "  absent   $label"
  fi
}
check /home/g/pepmlm_pipeline/data/targets/NCAM1.fasta 1 "NCAM1 target fasta"
check /home/g/pepmlm_pipeline/models/PepMLM-650M/pytorch_model.bin 1 "PepMLM weights"
check /home/g/pepmlm_pipeline/models/PepMLM-650M/config.json 1 "PepMLM config"
check /home/g/pepmlm_pipeline/models/PepMLM-650M/vocab.txt 1 "PepMLM vocab"
check /home/g/pepmlm_pipeline/models/esm2_t33_650M_UR50D/pytorch_model.bin 1 "ESM2 HF-format weights"
check /home/g/pepmlm_pipeline/models/moPPIt/modules/bindevaluator_modules/models.py 1 "moPPIt modules"
check /home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt 0 "BindEvaluator ckpt"
check /home/g/pepmlm_pipeline/models/esmfold_v1/pytorch_model.bin 0 "ESMFold weights"
check /home/g/GraphPep_v1.3/ckpt/model_param.ckpt 1 "GraphPep ckpt"
check /home/g/GraphPep_v1.3/bin/predict.py 1 "GraphPep predict.py"
echo
echo "exit_flag=$some_fail"
