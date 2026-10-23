#!/usr/bin/env bash
U="https://hf-mirror.com/ChatterjeeLab/PepMLM-650M/resolve/main/pytorch_model.bin"
echo "===== 1. headers (no -L) ====="
timeout 45 curl -4 -sS -D - -o /dev/null "$U" 2>&1 | head -20
echo
echo "===== 2. follow redirect, 1MB range, time it ====="
timeout 120 curl -4 -sSL -o /tmp/pep_head.bin -w "code=%{http_code} bytes=%{size_download} time=%{time_total} speed=%{speed_download}\n" -r 0-1048575 "$U" 2>&1
ls -la /tmp/pep_head.bin 2>&1
echo
echo "===== 3. same for a mid-file range (tests sustained transfer) ====="
timeout 120 curl -4 -sSL -o /tmp/pep_mid.bin -w "code=%{http_code} bytes=%{size_download} time=%{time_total} speed=%{speed_download}\n" -r 1000000000-1001048575 "$U" 2>&1
echo
echo "===== 4. HF_HUB direct (bypass mirror) reachable? ====="
timeout 30 curl -4 -sS -o /dev/null -w "hf.co direct code=%{http_code}\n" "https://huggingface.co/ChatterjeeLab/PepMLM-650M/resolve/main/config.json" 2>&1
echo
echo "===== 5. modelscope mirror? ====="
timeout 30 curl -4 -sS -o /dev/null -w "modelscope code=%{http_code}\n" "https://www.modelscope.cn" 2>&1
