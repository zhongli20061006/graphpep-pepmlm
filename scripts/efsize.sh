#!/usr/bin/env bash
M="https://hf-mirror.com"
U="$M/facebook/esmfold_v1/resolve/main/pytorch_model.bin"
echo "=== ESMFold pytorch_model.bin size (via redirect) ==="
# hf-mirror 302-redirects to the CAS backend; read Content-Length from the final response
curl -4 -sSL -o /dev/null -D /tmp/efh.txt --max-time 60 "$U" 2>&1
grep -iE 'content-length|x-linked-size|content-range' /tmp/efh.txt | tail -5
echo
echo "=== HEAD only (no body) ==="
curl -4 -sS -I -L --max-time 60 "$U" 2>&1 | grep -iE '^HTTP|content-length|x-linked-size' | head -10
echo
echo "=== BindEvaluator ckpt expected size ==="
curl -4 -sS -I -L --max-time 60 "$M/ChatterjeeLab/moPPIt/resolve/main/classifier_ckpt/finetuned_BindEvaluator.ckpt" 2>&1 | grep -iE '^HTTP|content-length|x-linked-size' | head -10
echo
echo "=== current local sizes ==="
du -sh /home/g/pepmlm_pipeline/models/esmfold_v1 /home/g/pepmlm_pipeline/models/moPPIt
stat -c '%n %s' /home/g/pepmlm_pipeline/models/esmfold_v1/pytorch_model.bin.part \
                /home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt.part 2>/dev/null
