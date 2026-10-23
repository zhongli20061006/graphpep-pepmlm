#!/usr/bin/env bash
export HF_ENDPOINT=https://hf-mirror.com
echo "===== A. hf-mirror: esm2 650M file HEAD ====="
timeout 40 curl -4 -sSL -o /dev/null -w "http=%{http_code} size=%{size_download} len=%{size_header}\n" -r 0-100 "https://hf-mirror.com/facebook/esm2_t33_650M_UR50D/resolve/main/pytorch_model.bin" 2>&1
echo
echo "===== B. hf-mirror api search PepMLM ====="
timeout 40 curl -4 -sS "https://hf-mirror.com/api/models?search=PepMLM&limit=20" 2>&1 | head -c 2000
echo; echo
echo "===== C. hf-mirror api search moPPIt ====="
timeout 40 curl -4 -sS "https://hf-mirror.com/api/models?search=moPPIt&limit=20" 2>&1 | head -c 2000
echo; echo
echo "===== D. hf-mirror api search moppit lowercase ====="
timeout 40 curl -4 -sS "https://hf-mirror.com/api/models?search=moppit&limit=20" 2>&1 | head -c 1000
echo; echo
echo "===== E. hf-mirror api search pepmlm lowercase ====="
timeout 40 curl -4 -sS "https://hf-mirror.com/api/models?search=pepmlm&limit=20" 2>&1 | head -c 1000
echo; echo
echo "===== F. huggingface_hub installed? ====="
/home/g/miniconda3/envs/graphpep_env/bin/pip list 2>/dev/null | grep -iE "huggingface|transformers|tokenizers|safetensors|esm" 2>&1
echo "(end)"
