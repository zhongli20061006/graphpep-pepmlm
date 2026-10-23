#!/usr/bin/env bash
export HF_ENDPOINT=https://hf-mirror.com
echo "=== PepMLM repo file list ==="
curl -4 -sSL "https://hf-mirror.com/api/models/ChatterjeeLab/PepMLM-650M" \
  | tr ',' '\n' | grep -o '"rfilename":"[^"]*"' | head -30
echo
echo "=== README ==="
curl -4 -sSL "https://hf-mirror.com/ChatterjeeLab/PepMLM-650M/resolve/main/README.md" | head -100
