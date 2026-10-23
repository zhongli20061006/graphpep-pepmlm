#!/usr/bin/env bash
echo "=== parallel dl log ==="
tail -14 /home/g/pepmlm_pipeline/logs/download_parallel.log 2>/dev/null
echo
echo "=== sizes ==="
du -sh /home/g/pepmlm_pipeline/models/* 2>/dev/null
echo
P=/home/g/pepmlm_pipeline/models/PepMLM-650M/pytorch_model.bin
B=/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt.part
E=/home/g/pepmlm_pipeline/models/esmfold_v1/pytorch_model.bin.part
size() { [ -f "$1" ] && stat -c %s "$1" || echo 0; }
p1=$(size "$P"); b1=$(size "$B"); e1=$(size "$E")
sleep 12
p2=$(size "$P"); b2=$(size "$B"); e2=$(size "$E")
echo "throughput over 12s:"
echo "  PepMLM      $(( (p2-p1)/12/1024 )) KB/s   total=$(( p2/1048576 )) MB"
echo "  BindEval    $(( (b2-b1)/12/1024 )) KB/s   total=$(( b2/1048576 )) MB"
echo "  ESMFold     $(( (e2-e1)/12/1024 )) KB/s   total=$(( e2/1048576 )) MB"
echo "  aggregate   $(( ((p2-p1)+(b2-b1)+(e2-e1))/12/1024 )) KB/s"
