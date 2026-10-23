#!/usr/bin/env bash
echo "===== A. curl -4 huggingface ====="
timeout 40 curl -4 -sS -o /dev/null -w "http_code=%{http_code} time=%{time_total}\n" https://huggingface.co 2>&1; echo "exit=$?"
echo
echo "===== B. curl -4 HF api (model meta) ====="
timeout 40 curl -4 -sS -w "\nhttp_code=%{http_code}\n" "https://huggingface.co/api/models/facebook/esm2_t33_650M_UR50D" 2>&1 | tail -5
echo
echo "===== C. search HF for PepMLM / moPPIt ====="
timeout 40 curl -4 -sS "https://huggingface.co/api/models?search=PepMLM&limit=10" 2>&1 | head -c 1500
echo
echo "--- moPPIt ---"
timeout 40 curl -4 -sS "https://huggingface.co/api/models?search=moPPIt&limit=10" 2>&1 | head -c 1500
echo
echo "===== D. npm-style: is hf-mirror reachable ====="
timeout 30 curl -4 -sS -o /dev/null -w "hf-mirror -> %{http_code}\n" https://hf-mirror.com 2>&1
echo
echo "===== E. github reachable (for esmfold source) ====="
timeout 30 curl -4 -sS -o /dev/null -w "github -> %{http_code}\n" https://github.com 2>&1
