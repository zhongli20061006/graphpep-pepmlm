#!/usr/bin/env bash
# 1) Measure per-connection vs aggregate throughput.
# 2) If per-connection is throttled, fetch the remaining BindEvaluator bytes with N
#    parallel range requests and CONCATENATE them ourselves.
# The existing .part is a valid prefix (curl -C - appended sequentially), so we resume
# from its length and then verify the assembled size exactly.
set -u
M="https://hf-mirror.com"
D=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt
FINAL=$D/finetuned_BindEvaluator.ckpt
PART=$FINAL.part
EXPECT=2832618887
PARTS=$D/ranges
mkdir -p "$PARTS"
URL="$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt"

echo "=== 1. throughput probe: 3 parallel 4MB ranges for 10s each ==="
for i in 0 1 2; do
  off=$(( 20000000 + i * 50000000 ))
  ( curl -4 -sS -o /dev/null -w "  stream$i: %{size_download} bytes in %{time_total}s = %{speed_download} B/s\n" \
      -r ${off}-$((off+4194303)) --max-time 20 "$URL" ) &
done
wait

echo
HAVE=0
[ -f "$PART" ] && HAVE=$(stat -c %s "$PART")
echo "=== 2. resume plan ==="
echo "  existing prefix : $HAVE bytes"
echo "  remaining       : $(( EXPECT - HAVE )) bytes"
if [ "$HAVE" -ge "$EXPECT" ]; then
  echo "  already complete; nothing to do"; exit 0
fi

# split remaining range across 6 workers
N=6
REM=$(( EXPECT - HAVE ))
CHUNK=$(( (REM + N - 1) / N ))
echo "  using $N parallel range requests of ~$CHUNK bytes"
pids=""
for i in $(seq 0 $((N-1))); do
  s=$(( HAVE + i * CHUNK ))
  e=$(( s + CHUNK - 1 ))
  [ "$e" -ge "$EXPECT" ] && e=$(( EXPECT - 1 ))
  [ "$s" -gt "$e" ] && continue
  out=$PARTS/part_$i
  echo "    worker $i: $s-$e"
  ( curl -4 -sSL --retry 10 --retry-delay 3 --retry-all-errors -C - -o "$out" \
      -r ${s}-${e} "$URL" ) &
  pids="$pids $!"
done
echo "  waiting for workers ..."
for p in $pids; do wait $p; done
echo "  workers done"
ls -la "$PARTS"

echo
echo "=== 3. assemble ==="
cp -f "$PART" "$FINAL.assembling" 2>/dev/null || true
for i in $(seq 0 $((N-1))); do
  out=$PARTS/part_$i
  [ -f "$out" ] && cat "$out" >> "$FINAL.assembling"
done
sz=$(stat -c %s "$FINAL.assembling" 2>/dev/null || echo 0)
echo "  assembled size: $sz (expected $EXPECT)"
if [ "$sz" -eq "$EXPECT" ]; then
  mv -f "$FINAL.assembling" "$FINAL"
  rm -rf "$PARTS" "$PART"
  echo "  SIZE OK -> $FINAL"
else
  echo "  SIZE MISMATCH -> keeping .assembling for inspection, NOT promoting"
fi
