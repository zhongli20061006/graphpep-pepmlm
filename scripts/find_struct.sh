#!/usr/bin/env bash
# Find experimental NCAM1 structures in the PDB (RCSB search API) and check for
# NCAM1-peptide complexes. Used for the Option-B receptor structure.
set -u
echo "=== RCSB search: NCAM1 human (P13591) X-ray/cryoEM structures ==="
curl -4 -sS --max-time 90 -X POST "https://search.rcsb.org/rcsbsearch/v2/query" \
  -H "Content-Type: application/json" \
  -d '{
    "query": {
      "type": "terminal",
      "service": "text",
      "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator": "exact_match",
        "value": "P13591"
      }
    },
    "return_type": "entry",
    "request_options": {"paginate": {"start": 0, "rows": 100},
                        "results_content_type": ["experimental"],
                        "sort": [{"sort_by": "rcsb_entry_info.resolution_combined",
                                  "direction": "asc"}]}
  }' 2>&1 | head -c 3000
echo
echo
echo "=== simpler: UniProt->PDB mapping from UniProt API ==="
curl -4 -sS --max-time 90 \
  "https://rest.uniprot.org/uniprotkb/P13591.json?fields=xref_pdb" 2>&1 \
  | python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
except Exception as e:
    print('parse fail', e); sys.exit()
pdb = []
for x in d.get('uniProtKBCrossReferences', []):
    if x.get('database') == 'PDB':
        pdb.append(x.get('id'))
print('n PDB entries:', len(pdb))
print(pdb)
" 2>&1
