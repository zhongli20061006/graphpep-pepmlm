#!/usr/bin/env bash
# Diagnose why the chunked download got 0 bytes.
set -u
M="https://hf-mirror.com"
URL="$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt"

echo "=== A. plain 1 MB GET (no -C -, no -r) ==="
timeout 60 curl -4 -sS -o /tmp/t1.bin -w "  code=%{http_code} bytes=%{size_download} time=%{time_total} speed=%{speed_download}\n" \
  -r 0-1048575 "$URL" 2>&1
echo "  file size: $(stat -c %s /tmp/t1.bin 2>/dev/null || echo none)"

echo
echo "=== B. same but through the redirect explicitly (-L) ==="
timeout 60 curl -4 -sSL -o /tmp/t2.bin -w "  code=%{http_code} bytes=%{size_download} time=%{time_total} speed=%{speed_download}\n" \
  -r 0-1048575 "$URL" 2>&1
echo "  file size: $(stat -c %s /tmp/t2.bin 2>/dev/null || echo none)"

echo
echo "=== C. the -C - plus -r combination that the script used (reproduce the bug) ==="
rm -f /tmp/t3.bin
# first, create a 1 MB partial
timeout 60 curl -4 -sSL -o /tmp/t3.bin -r 0-1048575 "$URL" 2>&1
echo "  after first fetch: $(stat -c %s /tmp/t3.bin 2>/dev/null || echo none) bytes"
# now resume: -C - with the SAME absolute -r window
timeout 60 curl -4 -sSL -C - -r 0-2097151 -o /tmp/t3.bin "$URL" 2>&1 | head -3
echo "  after -C - -r 0-2097151: $(stat -c %s /tmp/t3.bin 2>/dev/null || echo none) bytes"
echo "  (curl adds the existing size to the -r start, so it asked for 1048576-3145727)"

echo
echo "=== D. headers of the resource ==="
timeout 45 curl -4 -sSI -L "$URL" 2>&1 | grep -iE '^HTTP|content-length|accept-ranges|x-linked-size' | head -8

echo
echo "=== E. does the xet/older direct host still work? ==="
timeout 45 curl -4 -sS -o /dev/null -w "  hf-mirror root: %{http_code}\n" https://hf-mirror.com 2>&1

echo
echo "=== F. current DNS/network sanity ==="
timeout 20 getent hosts hf-mirror.com | head -3
