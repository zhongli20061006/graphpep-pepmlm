#!/usr/bin/env bash
# STOP all competing downloads and report the true state.
# Rationale: two ranged_dl_bindeval.sh instances plus a legacy hf_hub downloader were all
# writing the same target, so no single writer could make reliable progress.
set -u
D=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt

echo "=== killing ranged_dl + curl writers ==="
for pid in $(ps -ef | awk '/[r]anged_dl_bindeval|[c]url/ {print $2}'); do
  kill -9 "$pid" 2>/dev/null && echo "  killed $pid"
done
sleep 3

echo
echo "=== remaining writers ==="
ps -ef | grep -E '[r]anged_dl|[c]url|[d]ownload_models' || echo "  none"

echo
echo "=== disk state ==="
du -sh /home/g/pepmlm_pipeline/models/moPPIt 2>/dev/null
echo "ranges dir: $([ -d "$D/ranges" ] && du -sh "$D/ranges" || echo absent)"
echo "final ckpt: $([ -f "$D/finetuned_BindEvaluator.ckpt" ] && stat -c %s "$D/finetuned_BindEvaluator.ckpt" || echo MISSING)"
find /home/g/pepmlm_pipeline/models/moPPIt -name '*.incomplete' 2>/dev/null | while read f; do
  echo "incomplete: $(stat -c %s "$f")  $f"
done
