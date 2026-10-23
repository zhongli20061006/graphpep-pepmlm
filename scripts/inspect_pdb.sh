#!/usr/bin/env bash
# Inspect the 7 NCAM1 PDB entries: chains, lengths, ligands, peptides, resolution.
set -u
mkdir -p /home/g/pepmlm_pipeline/data/pdb
cd /home/g/pepmlm_pipeline/data/pdb
for id in 2E3V 2HAZ 2VKW 2VKX 3MTR 5AEA 5LKN; do
  low=$(echo "$id" | tr 'A-Z' 'a-z')
  [ -s "$low.pdb" ] || curl -4 -sSL --max-time 120 -o "$low.pdb" "https://files.rcsb.org/download/$id.pdb"
  echo "=============================================================="
  echo "=== $id  ($(stat -c %s "$low.pdb") bytes) ==="
  # title + resolution
  grep -m2 "^TITLE" "$low.pdb" | cut -c1-100
  grep -m1 "REMARK   2 RESOLUTION" "$low.pdb" | cut -c1-80
  # chain summary from SEQRES
  awk '$1=="SEQRES"{ch[$3]++} END{for(c in ch) printf "  chain %s: %d residues\n", c, ch[c]}' "$low.pdb" | sort
  # HETATM ligands (non-water)
  echo -n "  hetero (non-water, top): "
  awk '$1=="HETATM" && $4!="HOH"{print $4}' "$low.pdb" | sort -u | tr '\n' ' '
  echo
  # short chains = candidate peptides (<40 residues)
  awk '$1=="SEQRES"{n[$3]++} END{for(c in n) if(n[c]<40) printf "  SHORT CHAIN %s: %d residues (peptide?)\n", c, n[c]}' "$low.pdb"
done
