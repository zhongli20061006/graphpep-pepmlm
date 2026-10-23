#!/usr/bin/env bash
D=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt
echo "=== chunked download status ($(date +%H:%M:%S)) ==="
tail -6 /home/g/pepmlm_pipeline/logs/chunked_dl.log
echo
total=0
for f in "$D"/chunks/c_*; do [ -f "$f" ] && total=$(( total + $(stat -c %s "$f") )); done
echo "bytes on disk : $total / 2832618887  ($(( 100 * total / 2832618887 ))%)"
echo "chunks present: $(ls "$D"/chunks 2>/dev/null | wc -l) / 113"
echo "curl running  : $(ps -ef | grep -c '[c]url')"
echo "final present : $([ -f "$D/finetuned_BindEvaluator.ckpt" ] && stat -c %s "$D/finetuned_BindEvaluator.ckpt" || echo NO)"
