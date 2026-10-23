#!/usr/bin/env bash
# Are any of the NCAM1 PDB entries actually protein-protein COMPLEXES?
# If so, a bound peptide/complex partner gives a real interface to place candidates at.
set -u
cd /home/g/pepmlm_pipeline/data/pdb
for f in 2e3v.pdb 2haz.pdb 2vkw.pdb 2vkx.pdb 3mtr.pdb 5aea.pdb 5lkn.pdb; do
  [ -s "$f" ] || continue
  echo "=============================================================="
  echo "=== $f ==="
  # distinct polymer chains + their residue counts (from ATOM records)
  awk '$1=="ATOM"{
        ch=substr($0,22,1);
        key=ch":"substr($0,23,5);
        if(!(key in seen)){seen[key]=1; cnt[ch]++}
      }
      END{for(c in cnt) printf "  chain %s: %d residues\n", c, cnt[c]}' "$f" | sort
  # COMPND descriptions
  grep "^COMPND" "$f" | cut -c1-95 | head -6
done
