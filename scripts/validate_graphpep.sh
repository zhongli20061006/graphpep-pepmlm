#!/usr/bin/env bash
# Validate GraphPep end-to-end on its own bundled example (proves the invocation + outputs).
set -u
GP=/home/g/GraphPep_v1.3
echo "=== GraphPep example contents ==="
ls -la $GP/example/ 2>&1
find $GP/example -maxdepth 2 -type f | head -20
echo
echo "=== bin/ contents ==="
ls $GP/bin/
echo
echo "=== pick a case and run GraphPep with -each (-outeach) ==="
CASE=$(find $GP/example -name 'protein.pdb' | head -1)
if [ -z "$CASE" ]; then echo "no example protein.pdb found"; exit 1; fi
CDIR=$(dirname "$CASE")
echo "case dir: $CDIR"
ls -la "$CDIR"
DEC=$(find "$CDIR" -name '*decoy*' | head -1)
echo "decoys: $DEC"
echo
WORK=/tmp/gp_validate
rm -rf $WORK && mkdir -p $WORK
cd $CDIR
timeout 1800 bash $GP/GraphPep.sh \
    "$(realpath $(basename protein.pdb))" "$(realpath $(basename $DEC))" \
    -out $WORK/score.csv -each -outeach $WORK/score_each.csv 2>&1 | tail -25
echo
echo "=== outputs ==="
ls -la $WORK/
echo
echo "--- score.csv (head) ---"
head -6 $WORK/score.csv 2>&1
echo "--- score_each.csv (head) ---"
head -8 $WORK/score_each.csv 2>&1
