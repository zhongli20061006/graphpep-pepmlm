#!/usr/bin/env bash
# Retry the failed PDB downloads (connection resets) and parse them properly.
set -u
cd /home/g/pepmlm_pipeline/data/pdb
for id in 2E3V 5AEA 5LKN; do
  low=$(echo "$id" | tr 'A-Z' 'a-z')
  if [ -s "$low.pdb" ]; then echo "$id already present"; continue; fi
  for try in 1 2 3 4 5; do
    curl -4 -sSL --max-time 180 --retry 3 --retry-all-errors -o "$low.pdb" \
      "https://files.rcsb.org/download/$id.pdb" 2>/dev/null
    if [ -s "$low.pdb" ]; then break; fi
    echo "  $id attempt $try failed; retrying"
    sleep 3
  done
  echo "$id -> $(stat -c %s "$low.pdb" 2>/dev/null || echo MISSING) bytes"
done
echo
echo "=== summary of ALL entries ==="
for f in *.pdb; do
  [ -s "$f" ] || continue
  id=$(echo "$f" | tr 'a-z' 'A-Z' | sed 's/.PDB//')
  title=$(grep -m1 "^TITLE" "$f" | cut -c11-90)
  res=$(grep -m1 "REMARK   2 RESOLUTION" "$f" | awk '{print $4}')
  echo "--- $id  res=${res:-?}A"
  echo "    ${title:-<no title>}"
  # residue counts per chain from ATOM records (more reliable than SEQRES counts)
  awk '$1=="ATOM" && substr($0,22,1)!=" "{print substr($0,22,1)}' "$f" | sort | uniq -c \
    | awk '{printf "    chain %s: %s atoms\n", $2, $1}'
done
