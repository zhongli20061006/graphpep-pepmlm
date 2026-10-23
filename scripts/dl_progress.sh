#!/usr/bin/env bash
D=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt
echo "--- parts ---"
ls -la "$D/ranges" 2>/dev/null | tail -16
total=0
for f in "$D"/ranges/part_*; do
  [ -f "$f" ] || continue
  total=$(( total + $(stat -c %s "$f") ))
done
echo
echo "total so far: $total bytes = $(( total / 1048576 )) MB of 2701 MB"
echo "final present: $([ -f "$D/finetuned_BindEvaluator.ckpt" ] && stat -c %s "$D/finetuned_BindEvaluator.ckpt" || echo NO)"
echo "curl procs: $(ps -ef | grep -c '[c]url')"
