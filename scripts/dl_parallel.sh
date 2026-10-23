#!/usr/bin/env bash
# Parallel downloader for the remaining artifacts.
# Safe to run alongside dl3.sh: skips any file that is already fully present.
set -u
ROOT=/home/g/pepmlm_pipeline
MODELS=$ROOT/models
LOG=$ROOT/logs/download_parallel.log

M="https://hf-mirror.com"
BE_FINAL="$MODELS/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt"
BE_TMP="$MODELS/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt.part"
EF_FINAL="$MODELS/esmfold_v1/pytorch_model.bin"
EF_TMP="$MODELS/esmfold_v1/pytorch_model.bin.part"
EF_CFG="$MODELS/esmfold_v1/config.json"

mkdir -p "$(dirname "$BE_FINAL")" "$(dirname "$EF_FINAL")"

log() { echo "[$(date +%H:%M:%S)] $*"; }

fetch() {  # fetch <url> <final> <tmp> <minbytes> <label>
  local url="$1" final="$2" tmp="$3" minb="$4" label="$5"
  local have=0
  [ -f "$final" ] && have=$(stat -c '%s' "$final")
  if [ "$have" -ge "$minb" ]; then
    log "SKIP $label (already ${have} bytes)"
    return 0
  fi
  if [ -f "$tmp" ]; then
    log "RESUME $label from $(stat -c '%s' "$tmp") bytes"
  else
    log "START $label"
  fi
  timeout 10800 curl -4 -sSL --retry 8 --retry-delay 5 --retry-all-errors \
      -C - -o "$tmp" "$url" -w "  $label: code=%{http_code} bytes=%{size_download} time=%{time_total}s\n"
  local rc=$?
  local sz=0
  [ -f "$tmp" ] && sz=$(stat -c '%s' "$tmp")
  log "END $label rc=$rc size=$sz"
  if [ "$sz" -ge "$minb" ]; then
    mv -f "$tmp" "$final"
    log "COMMITTED $label -> $final"
  else
    log "INCOMPLETE $label ($sz < $minb); keeping .part for resume"
  fi
}

# --- BindEvaluator checkpoint (~2.83 GB) ---
(
  fetch "$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt" \
        "$BE_FINAL" "$BE_TMP" 2800000000 "BindEvaluator"
) >> "$LOG" 2>&1 &
BE_PID=$!

# --- ESMFold weights (~3 GB, exact size checked after) ---
(
  fetch "$M/facebook/esmfold_v1/resolve/main/config.json" "$EF_CFG" "$EF_CFG.part" 100 "esmfold_cfg"
  fetch "$M/facebook/esmfold_v1/resolve/main/pytorch_model.bin" "$EF_FINAL" "$EF_TMP" 1000000000 "esmfold_weights"
) >> "$LOG" 2>&1 &
EF_PID=$!

log "launched parallel fetchers: BindEvaluator pid=$BE_PID esmfold pid=$EF_PID"
wait $BE_PID $EF_PID
log "all parallel fetches finished"
du -sh "$MODELS"/* 2>&1
