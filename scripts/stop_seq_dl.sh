#!/usr/bin/env bash
# Stop the sequential BindEvaluator downloader (single writer on the .part file),
# then hand over to the multi-connection fetcher.
kill 2880 2878 2873 2869 2>/dev/null
sleep 3
echo "remaining curl procs: $(ps -ef | grep -c '[c]url')"
echo "part size now: $(stat -c %s /home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt/finetuned_BindEvaluator.ckpt.part 2>/dev/null || echo 0)"
