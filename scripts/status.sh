#!/usr/bin/env bash
echo "=== STEP1 progress ==="
echo "ppl cache entries:"
/home/g/miniconda3/envs/graphpep_env/bin/python - <<'PY'
import json, os
p = "/home/g/pepmlm_pipeline/cache/ppl_cache.json"
if os.path.isfile(p):
    d = json.load(open(p))
    print("  ", len(d), "peptides with cached PPL")
else:
    print("   (no cache yet)")
PY
echo
echo "=== output dir ==="
ls -la /home/g/pepmlm_pipeline/output/ 2>/dev/null
echo
echo "=== BindEvaluator download ==="
F=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt.part
G=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt
[ -f "$F" ] && echo "  .part = $(stat -c %s "$F") bytes ($(( $(stat -c %s "$F") / 1048576 )) MB of 2701 MB)"
[ -f "$G" ] && echo "  FINAL = $(stat -c %s "$G") bytes"
echo
echo "=== GPU ==="
nvidia-smi --query-gpu=memory.used,utilization.gpu --format=csv,noheader
echo
echo "=== running python ==="
ps -ef | grep -E 'step1|step2|prof' | grep -v grep
