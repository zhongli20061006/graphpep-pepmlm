#!/usr/bin/env bash
# Validate GraphPep end-to-end (corrected: inputs must be RELATIVE -- GraphPep symlinks them).
set -u
GP=/home/g/GraphPep_v1.3
CDIR=$GP/example/1B9J
WORK=/tmp/gp_validate
rm -rf $WORK && mkdir -p $WORK

echo "=== existing reference outputs in the repo ==="
echo "--- 1B9J/score.csv ---"; cat $CDIR/score.csv
echo "--- 1B9J/score_each.csv (head) ---"; head -6 $CDIR/score_each.csv
echo "--- 1B9J/decoys.pdb (MODEL count) ---"; grep -c "^MODEL" $CDIR/decoys.pdb
echo

cd "$CDIR"
echo "=== running GraphPep (relative inputs) ==="
timeout 1800 bash $GP/GraphPep.sh protein.pdb decoys.pdb \
    -out $WORK/score.csv -each -outeach $WORK/score_each.csv 2>&1 | tail -15
echo
echo "=== outputs ==="
ls -la $WORK/
echo
echo "--- fresh score.csv ---"
cat $WORK/score.csv 2>&1
echo
echo "--- fresh score_each.csv (head) ---"
head -8 $WORK/score_each.csv 2>&1
echo
echo "=== compare against the repo's stored score.csv ==="
if [ -f "$WORK/score.csv" ]; then
  echo "repo:  $(tr '\n' ' ' < $CDIR/score.csv)"
  echo "fresh: $(tr '\n' ' ' < $WORK/score.csv)"
fi
