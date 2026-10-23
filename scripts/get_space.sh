#!/usr/bin/env bash
U1="https://hf-mirror.com/spaces/TianlaiChen/PepMLM/raw/main/app.py"
U2="https://hf-mirror.com/spaces/TianlaiChen/PepMLM/resolve/main/app.py"
for u in "$U1" "$U2"; do
  echo "=== trying $u ==="
  curl -4 -sSL "$u" -o /tmp/app.py -w "code=%{http_code} bytes=%{size_download}\n" || true
  if [ -s /tmp/app.py ]; then
    echo "--- got it ---"
    cat /tmp/app.py
    exit 0
  fi
done
echo "=== both failed; trying HF api resolve-cache ==="
curl -4 -sSL "https://hf-mirror.com/api/spaces/TianlaiChen/PepMLM" 2>&1 | head -c 800
