#!/usr/bin/env bash
# Robust chunked download of the BindEvaluator checkpoint.
#
# Lessons baked in:
#   * `-C -` is MUTUALLY EXCLUSIVE with `-r` in curl ("--continue-at is mutually exclusive
#     with --range"). The earlier design used both, so every request failed instantly.
#   * The mirror drops large transfers mid-flight, so chunks must be small AND resumable.
#
# Strategy: for each fixed 32 MB window [s,e], track how many bytes are already on disk and
# request ONLY the missing sub-window, appending with >> . A drop therefore costs only the
# bytes in flight, and no curl resume flag is ever needed.
set -u
M="https://hf-mirror.com"
D=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt
FINAL=$D/finetuned_BindEvaluator.ckpt
CHUNKS=$D/chunks
URL="$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt"
EXPECT=2832618887
CHUNK_MB=${1:-32}
JOBS=${2:-6}
CHUNK=$(( CHUNK_MB * 1024 * 1024 ))
NCHUNK=$(( (EXPECT + CHUNK - 1) / CHUNK ))

rm -rf "$CHUNKS"; mkdir -p "$CHUNKS"
rm -f "$D"/finetuned_BindEvaluator.ckpt.part "$D"/finetuned_BindEvaluator.ckpt.assembling

echo "=== $EXPECT bytes in $NCHUNK chunks of ${CHUNK_MB} MB, $JOBS parallel ==="
echo "start $(date -Is)"

fetch_chunk() {
  local i="$1"
  local s=$(( i * CHUNK ))
  local e=$(( s + CHUNK - 1 ))
  [ "$e" -ge "$EXPECT" ] && e=$(( EXPECT - 1 ))
  local want=$(( e - s + 1 ))
  local out="$CHUNKS/c_$(printf '%05d' $i)"
  local attempt=0
  while :; do
    local have=0
    [ -f "$out" ] && have=$(stat -c %s "$out")
    if [ "$have" -eq "$want" ]; then return 0; fi
    if [ "$have" -gt "$want" ]; then
      echo "  chunk $i overshoot ($have > $want) -> truncating"
      truncate -s "$want" "$out"
      return 0
    fi
    attempt=$(( attempt + 1 ))
    if [ "$attempt" -gt 60 ]; then
      echo "  chunk $i GAVE UP after $attempt attempts ($have/$want)"
      return 1
    fi
    # request only the MISSING bytes, append them
    local rs=$(( s + have ))
    curl -4 -sSL --max-time 300 --retry 2 --retry-delay 2 --retry-all-errors \
         -r "${rs}-${e}" "$URL" >> "$out" 2>/dev/null
    sleep 0.3
  done
}

i=0
while [ "$i" -lt "$NCHUNK" ]; do
  batch=""
  for _ in $(seq 1 "$JOBS"); do
    [ "$i" -ge "$NCHUNK" ] && break
    fetch_chunk "$i" &
    batch="$batch $!"
    i=$(( i + 1 ))
  done
  for p in $batch; do wait "$p" || true; done
  total=0
  for f in "$CHUNKS"/c_*; do [ -f "$f" ] && total=$(( total + $(stat -c %s "$f") )); done
  echo "  progress: $total/$EXPECT bytes ($(( 100 * total / EXPECT ))%)  $(date +%H:%M:%S)"
done

echo
echo "=== verifying chunk sizes ==="
total=0; bad=0
for i in $(seq 0 $((NCHUNK-1))); do
  f="$CHUNKS/c_$(printf '%05d' $i)"
  s=$(( i * CHUNK )); e=$(( s + CHUNK - 1 )); [ "$e" -ge "$EXPECT" ] && e=$(( EXPECT - 1 ))
  want=$(( e - s + 1 ))
  have=$([ -f "$f" ] && stat -c %s "$f" || echo 0)
  total=$(( total + have ))
  if [ "$have" -ne "$want" ]; then echo "  BAD chunk $i: $have != $want"; bad=$(( bad + 1 )); fi
done
echo "  total = $total (expected $EXPECT)  bad_chunks=$bad"
if [ "$total" -ne "$EXPECT" ] || [ "$bad" -ne 0 ]; then
  echo "ABORT: chunk sizes wrong"; exit 2
fi

echo
echo "=== assembling ==="
cat $(ls "$CHUNKS"/c_* | sort) > "$FINAL.assembling"
sz=$(stat -c %s "$FINAL.assembling")
echo "assembled = $sz"
[ "$sz" -ne "$EXPECT" ] && { echo "ABORT: size mismatch"; exit 3; }

echo
echo "=== validating with torch.load ==="
if /home/g/miniconda3/envs/graphpep_env/bin/python - "$FINAL.assembling" <<'PY'
import sys, torch
ck = torch.load(sys.argv[1], map_location="cpu", weights_only=False)
print("  loaded type:", type(ck).__name__)
if isinstance(ck, dict):
    print("  keys:", list(ck.keys()))
    sd = ck.get("state_dict", {})
    print("  state_dict tensors:", len(sd))
    for k in list(sd)[:8]:
        print("   ", k, tuple(sd[k].shape))
PY
then
  mv -f "$FINAL.assembling" "$FINAL"
  rm -rf "$CHUNKS"
  echo "PROMOTED -> $FINAL ($(stat -c %s "$FINAL") bytes)"
else
  echo "ABORT: torch.load failed; keeping .assembling"
  exit 4
fi
echo "done $(date -Is)"
