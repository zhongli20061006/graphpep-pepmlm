#!/usr/bin/env bash
# Robust downloader using curl (bypasses the hanging Xet backend).
# huggingface.co is unreachable here; hf-mirror.com works over IPv4.
set -u
ROOT=/home/g/pepmlm_pipeline
MODELS=$ROOT/models
mkdir -p "$MODELS" "$ROOT/logs"
cd "$MODELS"

M="https://hf-mirror.com"
get() {  # get <url> <outfile> <label>
  local url="$1" out="$2" label="$3"
  echo "[$(date +%H:%M:%S)] START $label -> $out"
  timeout 7200 curl -4 -sSL --retry 5 --retry-delay 5 --retry-all-errors \
      -C - -o "$out" -w "  $label: code=%{http_code} bytes=%{size_download} time=%{time_total}s speed=%{speed_download}B/s\n" \
      "$url" 2>&1
  echo "[$(date +%H:%M:%S)] DONE  $label size=$(du -h "$out" 2>/dev/null | cut -f1)"
}

echo "=================== downloads start $(date -Is) ==================="

# ---- 1. PepMLM-650M (transformers format) ----
mkdir -p PepMLM-650M
get "$M/ChatterjeeLab/PepMLM-650M/resolve/main/config.json"           PepMLM-650M/config.json
get "$M/ChatterjeeLab/PepMLM-650M/resolve/main/vocab.txt"             PepMLM-650M/vocab.txt
get "$M/ChatterjeeLab/PepMLM-650M/resolve/main/tokenizer_config.json" PepMLM-650M/tokenizer_config.json
get "$M/ChatterjeeLab/PepMLM-650M/resolve/main/special_tokens_map.json" PepMLM-650M/special_tokens_map.json
get "$M/ChatterjeeLab/PepMLM-650M/resolve/main/pytorch_model.bin"     PepMLM-650M/pytorch_model.bin

# ---- 2. BindEvaluator checkpoint (moPPIt) ----
mkdir -p moPPIt/classifier_ckpt
get "$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt" \
    moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt

echo "=================== PEPMLM+BINDEVALUATOR done $(date -Is) ==================="
du -sh PepMLM-650M moPPIt 2>&1
echo "=================== ESMFOLD next (large, ~3-14GB) ==================="

# ---- 3. ESMFold v1 ----
mkdir -p esmfold_v1
for f in config.json pytorch_model.bin; do
  get "$M/facebook/esmfold_v1/resolve/main/$f" "esmfold_v1/$f"
done

echo "=================== all downloads end $(date -Is) ==================="
du -sh "$MODELS"/* 2>&1
