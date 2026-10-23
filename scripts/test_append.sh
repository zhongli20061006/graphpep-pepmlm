#!/usr/bin/env bash
# Verify the append-resume fetch pattern on a small scale.
set -u
M="https://hf-mirror.com"
U="$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt"
cd /tmp
rm -f probe.bin
echo "=== two sequential 4MB append requests ==="
curl -4 -sSL --max-time 120 -r 0-4194303 "$U" >> probe.bin
echo "  after chunk 1: $(stat -c %s probe.bin) bytes (expect 4194304)"
curl -4 -sSL --max-time 120 -r 4194304-8388607 "$U" >> probe.bin
echo "  after chunk 2: $(stat -c %s probe.bin) bytes (expect 8388608)"
echo
echo "=== simulate a mid-chunk drop then resume the missing bytes ==="
rm -f probe2.bin
# fetch a deliberately truncated piece (first 1MB of the window)
curl -4 -sSL --max-time 120 -r 0-1048575 "$U" >> probe2.bin
echo "  partial: $(stat -c %s probe2.bin) bytes"
have=$(stat -c %s probe2.bin)
# resume from the exact byte after what we have, up to 4MB total
curl -4 -sSL --max-time 120 -r "${have}-4194303" "$U" >> probe2.bin
echo "  after resume: $(stat -c %s probe2.bin) bytes (expect 4194304)"
echo
echo "=== are the two files identical over the first 4MB? ==="
if cmp -s probe.bin <(head -c 4194304 probe2.bin); then
  echo "  IDENTICAL -> append-resume is byte-exact"
else
  echo "  DIFFER -> append-resume is WRONG"
fi
