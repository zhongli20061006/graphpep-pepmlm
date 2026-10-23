#!/usr/bin/env bash
# Definitively verify that append-resume reproduces the same bytes as a single clean fetch.
set -u
M="https://hf-mirror.com"
U="$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt"
cd /tmp

echo "=== reference: one clean 4 MB request ==="
rm -f ref.bin
curl -4 -sSL --max-time 180 -r 0-4194303 "$U" -o ref.bin
echo "  ref.bin = $(stat -c %s ref.bin) bytes"
md5sum ref.bin | awk '{print "  md5:", $1}'

echo
echo "=== split+append: 1 MB then resume the rest ==="
rm -f split.bin
curl -4 -sSL --max-time 180 -r 0-1048575 "$U" -o split.bin
echo "  first piece = $(stat -c %s split.bin) bytes"
have=$(stat -c %s split.bin)
curl -4 -sSL --max-time 180 -r "${have}-4194303" "$U" >> split.bin
echo "  after resume = $(stat -c %s split.bin) bytes"
md5sum split.bin | awk '{print "  md5:", $1}'

echo
echo "=== chunk-wise md5 comparison (1 MB blocks) ==="
ok=1
for blk in 0 1 2 3; do
  a=$(dd if=ref.bin   bs=1048576 skip=$blk count=1 2>/dev/null | md5sum | cut -d' ' -f1)
  b=$(dd if=split.bin bs=1048576 skip=$blk count=1 2>/dev/null | md5sum | cut -d' ' -f1)
  if [ "$a" = "$b" ]; then
    echo "  block $blk: MATCH  $a"
  else
    echo "  block $blk: DIFFER ref=$a split=$b"
    ok=0
  fi
done
echo
if [ "$ok" = "1" ]; then
  echo "RESULT: append-resume is BYTE-EXACT"
else
  echo "RESULT: append-resume is NOT byte-exact"
fi
