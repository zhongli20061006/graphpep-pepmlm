"""Bounded primary-source retrieval with local reuse; no sequence submission or models."""
from __future__ import annotations

import hashlib
import json
import urllib.request
from pathlib import Path

ROOT = Path('/mnt/e/Users/graphpep与pepmlm')
OUT = ROOT / 'output/structure_input_audit_bclw_20261020'
OLD = ROOT / 'output/bh3_bclw_source_audit_20261020/raw'
RAW = OUT / 'sources'
RAW.mkdir(exist_ok=True)
manifest: list[dict] = []


def fetch(url: str, filename: str, old: str | None = None, body: dict | None = None) -> bytes:
    """Reuse verified original material where available, else save a bounded request."""
    path = OLD / old if old else RAW / filename
    reused = bool(old)
    if not path.exists():
        req = urllib.request.Request(url, data=None if body is None else json.dumps(body).encode(),
                                     headers={'User-Agent': 'Bclw-input-audit/1.0', 'Content-Type': 'application/json'})
        with urllib.request.urlopen(req, timeout=45) as response:
            data = response.read(3_000_001)
        if len(data) > 3_000_000:
            raise ValueError('Per-source 3 MB audit bound exceeded')
        path.write_bytes(data)
    data = path.read_bytes()
    manifest.append({'url': url, 'path': str(path.relative_to(ROOT)), 'bytes': len(data),
                     'sha256': hashlib.sha256(data).hexdigest(), 'reused_prior_source': reused})
    return data


def main() -> None:
    query = {'query': {'type': 'terminal', 'service': 'text', 'parameters': {
        'attribute': 'rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession',
        'operator': 'exact_match', 'value': 'Q92843'}}, 'return_type': 'entry',
        'request_options': {'paginate': {'start': 0, 'rows': 30}}}
    data = json.loads(fetch('https://search.rcsb.org/rcsbsearch/v2/query', 'q92843_search.json', body=query))
    # 4K5A/B are explicit primary-paper references missed by Q92843 accession linkage.
    ids = sorted({r['identifier'].upper() for r in data.get('result_set', [])} | {'1ZY3', '4CIM', '3FDL', '4K5A', '4K5B'})
    assert len(ids) <= 30, 'Bounded audit needs manual review before expanding'
    summaries = []
    for pdb in ids:
        prefix = f'rcsb_{pdb.lower()}'
        old_entry = prefix + '_entry.json'
        entry = json.loads(fetch(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}', prefix + '_entry.json',
                                 old_entry if (OLD / old_entry).exists() else None))
        entities = []
        for eid in entry['rcsb_entry_container_identifiers'].get('polymer_entity_ids', []):
            filename = prefix + '_entity' + eid + '.json'
            entity = json.loads(fetch(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{eid}', filename,
                                      filename if (OLD / filename).exists() else None))
            ep = entity['entity_poly']
            entities.append({'entity_id': eid, 'description': entity['rcsb_polymer_entity'].get('pdbx_description'),
                             'chains': ep.get('pdbx_strand_id'), 'sequence': ep.get('pdbx_seq_one_letter_code_can', '').replace('\n', ''),
                             'length': ep.get('rcsb_sample_sequence_length'), 'mutation': entity['rcsb_polymer_entity'].get('pdbx_mutation'),
                             'references': entity['rcsb_polymer_entity_container_identifiers'].get('reference_sequence_identifiers'),
                             'align': entity.get('rcsb_polymer_entity_align', [])})
        summaries.append({'pdb_id': pdb, 'title': entry['struct']['title'], 'method': entry.get('exptl'),
                          'resolution': entry.get('rcsb_entry_info', {}).get('resolution_combined'),
                          'citation': entry.get('rcsb_primary_citation'), 'entities': entities})
    (OUT / 'structure_database_inventory.json').write_text(json.dumps({
        'query_scope': 'RCSB accession-linked Q92843 entries, at most 30, plus three handoff references; experimental only checked from exptl. Not a proof against unlinked/homologous structures.',
        'total_query_count': data.get('total_count'), 'entries': summaries}, ensure_ascii=False, indent=2), encoding='utf-8')
    for pdb in ('1ZY3', '4CIM', '3FDL'):
        fetch(f'https://files.rcsb.org/download/{pdb}.cif', pdb + '.cif')
    try:
        fetch('https://pmc.ncbi.nlm.nih.gov/articles/PMC4276867/?report=xml', 'Lee2014_4CIM_source.html')
    except Exception as error:
        manifest.append({'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4276867/?report=xml',
                         'status': 'failed_not_evidence', 'error': str(error)})
    manifest.append({'url': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4276867/fullTextXML',
                     'status': 'earlier_attempt_failed_not_evidence', 'error': 'HTTP 500; kept as failed retrieval; web tool primary article read separately'})
    (OUT / 'source_retrieval_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    for row in summaries:
        print(row['pdb_id'], row['title'])
        for e in row['entities']:
            print(' ', e['entity_id'], e['description'], e['chains'], e['length'], e['mutation'], e['sequence'])


if __name__ == '__main__':
    try:
        main()
    finally:
        (OUT / 'source_retrieval_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
