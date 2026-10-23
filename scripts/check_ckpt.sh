#!/usr/bin/env bash
D=/home/g/pepmlm_pipeline/models/moPPIt
echo "=== moPPIt tree sizes ==="
du -sh "$D"/* 2>/dev/null
echo
echo "=== BindEvaluator checkpoint candidates ==="
find "$D" -name '*.ckpt' -o -name '*.pt' -o -name '*.pth' 2>/dev/null | while read f; do
  printf "  %12s  %s\n" "$(stat -c %s "$f")" "$f"
done
echo
echo "=== expected sizes ==="
echo "  finetuned_BindEvaluator.ckpt expected = 2832618887"
echo
F="$D/classifier_ckpt/finetuned_BindEvaluator.ckpt"
if [ -f "$F" ]; then
  echo "FINAL EXISTS: $(stat -c %s "$F") bytes"
else
  echo "FINAL MISSING"
fi
echo
echo "=== ranged parts ==="
total=0
for f in "$D"/ranges/part_*; do
  [ -f "$f" ] || continue
  total=$(( total + $(stat -c %s "$f") ))
done
echo "  parts total = $total (expected 2832618887)"
echo "  curl procs   = $(ps -ef | grep -c '[c]url')"
echo
echo "=== partial downloads / incomplete ==="
find "$D" -name '*.incomplete' -o -name '*.part' 2>/dev/null | while read f; do
  printf "  %12s  %s\n" "$(stat -c %s "$f")" "$f"
done
