#!/usr/bin/env bash
# Fetch the whole BindEvaluator checkpoint as N clean byte ranges and concatenate.
#
# Why: the mirror throttles per connection (~200-500 KB/s) while aggregate bandwidth is
# much higher, so one sequential stream is the bottleneck. Ranges are requested from 0 so
# the assembly does not depend on any previously interrupted .part being a clean prefix.
#
# Safety: -r ONLY (no -C -: with -r, curl already knows the exact byte window, and -C -
# would resume from the output file's own size instead, producing a misaligned chunk).
# The assembled file must match the server's Content-Length exactly, and must load with
# torch.load, before it is promoted.
set -u
M="https://hf-mirror.com"
D=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt
FINAL=$D/finetuned_BindEvaluator.ckpt
PARTS=$D/ranges
URL="$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt"
EXPECT=2832618887
N=${1:-12}

rm -rf "$PARTS"; mkdir -p "$PARTS"
rm -f "$D"/finetuned_BindEvaluator.ckpt.part "$FINAL"

CHUNK=$(( (EXPECT + N - 1) / N ))
echo "=== fetching $EXPECT bytes in $N ranges of ~$CHUNK ==="
echo "start $(date -Is)"

pids=""
for i in $(seq 0 $((N-1))); do
  s=$(( i * CHUNK ))
  e=$(( s + CHUNK - 1 ))
  [ "$e" -ge "$EXPECT" ] && e=$(( EXPECT - 1 ))
  [ "$s" -gt "$e" ] && continue
  out=$PARTS/part_$(printf '%02d' $i)
  (
    for attempt in 1 2 3 4 5; do
      # -r ONLY; no -C -. Retry whole range on failure so the part is always complete.
      rm -f "$out"
      curl -4 -sSL --max-time 7200 -r "${s}-${e}" -o "$out" "$URL" \
        && [ "$(stat -c %s "$out")" -eq $(( e - s + 1 )) ] && exit 0
      echo "  range $i attempt $attempt incomplete, retrying"
      sleep 3
    done
    echo "  RANGE $i FAILED"
    exit 1
  ) &
  pids="$pids $!"
done

fail=0
for p in $pids; do wait "$p" || fail=1; done
echo "workers done $(date -Is) fail=$fail"

echo
echo "=== part sizes ==="
ls -la "$PARTS"
total=0
for f in "$PARTS"/part_*; do
  [ -f "$f" ] && total=$(( total + $(stat -c %s "$f") ))
done
echo "sum of parts = $total (expected $EXPECT)"

if [ "$total" -ne "$EXPECT" ]; then
  echo "ABORT: part sizes do not sum to the expected length; not assembling"
  exit 2
fi

echo
echo "=== assembling in order ==="
cat $(ls "$PARTS"/part_* | sort) > "$FINAL.assembling"
sz=$(stat -c %s "$FINAL.assembling")
echo "assembled = $sz"
if [ "$sz" -ne "$EXPECT" ]; then
  echo "ABORT: assembled size mismatch"
  exit 3
fi

echo
echo "=== verifying the checkpoint actually loads ==="
if /home/g/miniconda3/envs/graphpep_env/bin/python - "$FINAL.assembling" <<'PY'
import sys, torch
p = sys.argv[1]
ck = torch.load(p, map_location="cpu", weights_only=False)
print("  torch.load OK, type:", type(ck).__name__)
if isinstance(ck, dict):
    print("  top-level keys:", list(ck.keys())[:6])
    sd = ck.get("state_dict", ck)
    print("  n tensors:", len(sd))
    ks = list(sd.keys())
    print("  first keys:", ks[:4])
PY
then
  mv -f "$FINAL.assembling" "$FINAL"
  rm -rf "$PARTS"
  echo "PROMOTED -> $FINAL ($(stat -c %s "$FINAL") bytes)"
else
  echo "ABORT: checkpoint failed to load; leaving .assembling for inspection"
  exit 4
fi
echo "done $(date -Is)"
