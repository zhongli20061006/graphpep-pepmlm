#!/usr/bin/env bash
# Run GraphPep from a clean work dir, capturing ALL output and locating the outputs.
set -u
GP=/home/g/GraphPep_v1.3
WORK=/home/g/pepmlm_pipeline/work/gp_test        # avoid any 'tmp' in the path
rm -rf "$WORK" && mkdir -p "$WORK"
cd "$GP/example/1B9J"

echo "=== run: simple path, -out with a RELATIVE name (GraphPep copies to cwd) ==="
cp -f protein.pdb "$WORK/" && cp -f decoys.pdb "$WORK/"
cd "$WORK"
timeout 900 bash $GP/GraphPep.sh protein.pdb decoys.pdb -out score_simple.csv 2>&1 | tail -8
echo
echo "=== files in work dir ==="
ls -la "$WORK"
echo
echo "--- score_simple.csv ---"
cat "$WORK/score_simple.csv" 2>&1
echo
echo "=== reference ==="
cat "$GP/example/1B9J/score.csv"
echo
echo "=== now the -each path ==="
timeout 900 bash $GP/GraphPep.sh protein.pdb decoys.pdb -out score_e.csv -each -outeach score_each_fresh.csv 2>&1 | tail -12
echo
echo "--- files ---"; ls -la "$WORK"
echo "--- score_e.csv ---"; cat "$WORK/score_e.csv" 2>&1
echo "--- score_each_fresh.csv (head) ---"; head -6 "$WORK/score_each_fresh.csv" 2>&1
