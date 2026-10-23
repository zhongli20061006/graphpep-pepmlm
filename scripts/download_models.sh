#!/usr/bin/env bash
# Background downloader: fetch the three model families from the hf-mirror.com mirror.
# huggingface.co itself is unreachable from this WSL network (IPv6-only DNS, no v6 route).
set -u
export HF_ENDPOINT=https://hf-mirror.com
ROOT=/home/g/pepmlm_pipeline
MODELS=$ROOT/models
mkdir -p "$MODELS" "$ROOT/output" "$ROOT/data/targets" "$ROOT/scripts" "$ROOT/logs"
LOG=$ROOT/logs/download.log
exec > >(tee -a "$LOG") 2>&1

echo "=================== download start $(date -Is) ==================="

PY=/home/g/miniconda3/envs/graphpep_env/bin/python
PIP=/home/g/miniconda3/envs/graphpep_env/bin/pip

echo "--- [0/4] installing huggingface_hub + transformers ---"
timeout 900 $PIP install -q --disable-pip-version-check \
    "huggingface_hub[cli]" sentencepiece 2>&1 | tail -5
echo "hf_hub install exit: $?"
$PY -c "import huggingface_hub, sys; print('huggingface_hub', huggingface_hub.__version__)" 2>&1

echo
echo "--- [1/4] PepMLM-650M (2.6 GB) ---"
$PY - <<'PYEOF'
from huggingface_hub import snapshot_download
import os
os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")
p = snapshot_download(
    repo_id="ChatterjeeLab/PepMLM-650M",
    local_dir="/home/g/pepmlm_pipeline/models/PepMLM-650M",
    allow_patterns=["*.json", "*.txt", "*.bin", "*.safetensors", "*.model"],
    max_workers=4,
)
print("PepMLM ->", p)
PYEOF
echo "step1 exit: $?"

echo
echo "--- [2/4] moPPIt BindEvaluator (selective) ---"
$PY - <<'PYEOF'
from huggingface_hub import snapshot_download
p = snapshot_download(
    repo_id="ChatterjeeLab/moPPIt",
    local_dir="/home/g/pepmlm_pipeline/models/moPPIt",
    allow_patterns=[
        "bindevaluator.py", "moppit.py", "README.md",
        "scripts/*",
        "modules/**", "models/**", "utils/**", "flow_matching/**",
        "classifier_code/**",
        "classifier_ckpt/finetuned_BindEvaluator.ckpt",
        "classifier_ckpt/binding_affinity_pooled.pt",
        "classifier_ckpt/*.json",
        "ckpt/**",
    ],
    ignore_patterns=["*.pyc", "__pycache__/**"],
    max_workers=4,
)
print("moPPIt ->", p)
PYEOF
echo "step2 exit: $?"

echo
echo "--- [3/4] ESMFold v1 (esmfold_v1 pytorch_model.bin) ---"
$PY - <<'PYEOF'
from huggingface_hub import snapshot_download
try:
    p = snapshot_download(
        repo_id="facebook/esmfold_v1",
        local_dir="/home/g/pepmlm_pipeline/models/esmfold_v1",
        max_workers=4,
    )
    print("esmfold ->", p)
except Exception as e:
    print("ESMFOLD-DOWNLOAD-FAILED:", type(e).__name__, e)
PYEOF
echo "step3 exit: $?"

echo
echo "--- [4/4] verification ---"
du -sh "$MODELS"/* 2>&1
echo "=================== download end $(date -Is) ==================="
