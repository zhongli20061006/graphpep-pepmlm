"""Read existing pilot artifacts; compute audit geometry only, never inference/reward scores."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT = Path('/mnt/e/Users/graphpep与pepmlm')
OUT = ROOT / 'output/structure_input_audit_bclw_20261020'
PILOT = ROOT / 'output/bh3_bclw_six_complex_pilot_20261020'
SOURCE = ROOT / 'output/bh3_bclw_source_audit_20261020'


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def writej(name: str, value: object) -> None:
    (OUT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def writecsv(name: str, rows: list[dict]) -> None:
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def chains(path: Path) -> dict:
    model = PDBParser(QUIET=True).get_structure('audit', path)[0]
    result = {}
    for chain in model:
        residues = [r for r in chain if r.id[0] == ' ']
        result[chain.id] = {
            'sequence': ''.join(seq1(r.resname) for r in residues),
            'numbers': [r.id[1] for r in residues],
            'ca': np.array([r['CA'].coord for r in residues]),
            'plddt': np.array([r['CA'].bfactor for r in residues]),
            'atoms': np.array([a.coord for r in residues for a in r if a.element not in ('H', 'D')]),
            'atom_res': np.array([i for i, r in enumerate(residues) for a in r if a.element not in ('H', 'D')]),
        }
    return result


def geometry(chain_data: dict) -> dict:
    a, b = chain_data['A'], chain_data['B']
    dist2 = np.sum((a['atoms'][:, None, :] - b['atoms'][None, :, :]) ** 2, axis=-1)
    hit = np.argwhere(dist2 <= 36)
    pairs = sorted({(int(a['atom_res'][i]), int(b['atom_res'][j])) for i, j in hit})
    return {'minimum': float(np.sqrt(dist2.min())), 'clashes': int((dist2 < 2.25).sum()),
            'atom_pairs': len(hit), 'residue_pairs': pairs}


def main() -> None:
    original_manifest = readj(PILOT / 'run_result_manifest_v1.json')
    integrity = []
    for row in original_manifest['files']:
        path = PILOT / row['path']
        ok = path.exists() and path.stat().st_size == row['bytes'] and sha(path) == row['sha256']
        integrity.append({'path': row['path'], 'match': ok})
    writej('original_integrity_check.json', {'checked': len(integrity),
           'mismatch_count': sum(not r['match'] for r in integrity), 'rows': integrity})
    frozen = readj(PILOT / 'run_config_frozen_v1.json')
    records = {}
    header = None
    for line in (PILOT / 'complexes.fasta').read_text().splitlines():
        if line.startswith('>'):
            header = line[1:]
            records[header] = ''
        elif line.strip():
            records[header] += line.strip()
    uniprot = ''.join((SOURCE / 'raw/Q92843_BCLW.fasta').read_text().splitlines()[1:])[:164]
    expected = frozen['study_scope']['sequences']
    execution = readj(PILOT / 'run_execution_record_v1.json')
    config_rows, audit_rows, a3m_rows, log_rows = [], [], [], []
    old_gate = {(int(r['seed']), r['complex_id']): r for r in csv.DictReader(
        (PILOT / 'structure_gate_audit_v1.csv').open(encoding='utf-8-sig'))}
    raw_data = {}
    for seed in frozen['prediction']['seeds']:
        directory = PILOT / f'seed_{seed}/structures'
        config = readj(directory / 'config.json')
        event = next(e for e in execution['events'] if e['label'] == f'afm_seed_{seed}')
        checks = {'version': '1.6.3', 'msa_mode': 'single_sequence', 'model_type': 'alphafold2_multimer_v3',
                  'num_models': 1, 'num_recycles': 3, 'recycle_early_stop_tolerance': 0,
                  'num_ensemble': 1, 'model_order': [1], 'use_templates': False, 'num_relax': 0,
                  'pair_mode': 'unpaired', 'num_seeds': 1, 'random_seed': seed,
                  'use_fast_kernels': True, 'compile_mode': 'fast', 'max_seq': 1, 'max_extra_seq': 1}
        for key, wanted in checks.items():
            config_rows.append({'seed': seed, 'field': key, 'expected_effective': json.dumps(wanted),
                                'actual': json.dumps(config.get(key)), 'match': config.get(key) == wanted})
        log_text = (directory / 'log.txt').read_text()
        recycles = re.findall(r'(alphafold2_multimer_v3_model_1_seed_\d+) recycle=(\d+) pLDDT=([\d.]+) pTM=([\d.]+) ipTM=([\d.]+)', log_text)
        for i, (model, recycle, plddt, ptm, iptm) in enumerate(recycles):
            log_rows.append({'seed': seed, 'complex_id': list(records)[i // 4], 'recycle': int(recycle),
                             'whole_complex_plddt_log': float(plddt), 'ptm_log': float(ptm), 'iptm_log': float(iptm)})
        assert len(recycles) == 12 and event['returncode'] == 0
        for name, sequence in records.items():
            receptor, peptide = sequence.split(':')
            role = 'parent' if name.endswith('parent') else 'gain_candidate' if name.endswith('gain') else 'loss_candidate'
            assert receptor == uniprot and peptide == expected[role]
            a3m_path = directory / (name + '.a3m')
            a3m = a3m_path.read_text().splitlines()
            seqrows = [s for s in a3m if not s.startswith(('#', '>'))]
            # ColabFold serialized single-sequence unpaired MSA: one padded row per chain.
            a3m_ok = a3m[0] == '#164,23\t1,1' and seqrows == [receptor + '-' * 23, '-' * 164 + peptide]
            a3m_rows.append({'seed': seed, 'complex_id': name, 'header': a3m[0], 'serialized_rows': len(seqrows),
                             'receptor_biological_sequences': 1, 'peptide_biological_sequences': 1,
                             'homolog_sequences': 0, 'exact_query_and_gap_layout': a3m_ok, 'sha256': sha(a3m_path)})
            raw_path = next(directory.glob(name + '_unrelaxed_rank_001_*.pdb'))
            relaxed_path = next((PILOT / f'seed_{seed}/relaxed').glob(name + '_relaxed_rank_001_*.pdb'))
            score_path = next(directory.glob(name + '_scores_rank_001_*.json'))
            score = readj(score_path)
            raw, relaxed = chains(raw_path), chains(relaxed_path)
            raw_data[seed, name] = raw
            identity = all(d['A']['sequence'] == receptor and d['B']['sequence'] == peptide and
                           d['A']['numbers'] == list(range(1, 165)) and d['B']['numbers'] == list(range(1, 24))
                           for d in (raw, relaxed))
            plddt, pae = np.asarray(score['plddt']), np.asarray(score['pae'])
            assert plddt.shape == (187,) and pae.shape == (187, 187)
            pae_export = json.loads((directory / (name + '_predicted_aligned_error_v1.json')).read_text())
            if isinstance(pae_export, list):
                pae_export = pae_export[0]
            exported_matrix = np.asarray(pae_export['predicted_aligned_error'])
            export_difference = float(np.max(np.abs(exported_matrix - pae)))
            inter = np.concatenate((pae[:164, 164:].ravel(), pae[164:, :164].ravel()))
            rg, ag = geometry(raw), geometry(relaxed)
            gate = identity and score['iptm'] >= .7 and plddt[164:].mean() >= 70 and inter.max() <= 10 and ag['minimum'] >= 1.5 and len(ag['residue_pairs']) > 0
            prior = old_gate[seed, name]
            agreement = (abs(float(prior['iptm']) - score['iptm']) < 1e-9 and
                         abs(float(prior['peptide_mean_ca_plddt']) - float(plddt[164:].mean())) < 1e-8 and
                         abs(float(prior['interface_pae_max_A']) - float(inter.max())) < 1e-9 and
                         int(prior['unique_contact_residue_pairs_le_6A']) == len(ag['residue_pairs']) and
                         abs(float(prior['minimum_interchain_heavyatom_distance_A']) - ag['minimum']) < 1e-5 and
                         (prior['scoring_eligible'] == 'True') == bool(gate))
            contact_pae = [value for a, b in ag['residue_pairs'] for value in (pae[a, 164+b], pae[164+b, a])]
            row = {'seed': seed, 'complex_id': name, 'identity_raw_and_relaxed': identity,
                   'a3m_match': a3m_ok, 'iptm_json_rounded': score['iptm'], 'ptm_json_rounded': score['ptm'],
                   'receptor_mean_plddt': float(plddt[:164].mean()), 'receptor_fraction_plddt_ge70': float((plddt[:164] >= 70).mean()),
                   'peptide_mean_plddt': float(plddt[164:].mean()),
                   'pdb_ca_vs_json_max_abs_plddt_difference': float(np.max(np.abs(np.r_[raw['A']['plddt'], raw['B']['plddt']] - plddt))),
                   'exported_PAE_vs_scores_PAE_max_abs_difference': export_difference,
                   'all_directed_crosschain_pae_max_A': float(inter.max()),
                   'all_directed_crosschain_pae_median_A_DIAGNOSTIC_ONLY': float(np.median(inter)),
                   'contact_directed_pae_median_A_DIAGNOSTIC_ONLY': float(np.median(contact_pae)),
                   'receptor_intrachain_pae_median_A': float(np.median(pae[:164, :164])),
                   'raw_min_distance_A': rg['minimum'], 'relaxed_min_distance_A': ag['minimum'],
                   'relaxed_clash_count_lt1p5A': ag['clashes'], 'relaxed_residue_pairs_le6A': len(ag['residue_pairs']),
                   'relaxed_heavy_atom_pairs_le6A_DIAGNOSTIC_ONLY': ag['atom_pairs'],
                   'raw_to_relaxed_CA_RMSD_without_alignment_A': float(np.sqrt(np.mean(np.sum((np.r_[raw['A']['ca'], raw['B']['ca']] - np.r_[relaxed['A']['ca'], relaxed['B']['ca']]) ** 2, axis=-1)))),
                   'frozen_gate_pass': bool(gate), 'prior_gate_table_agreement': agreement,
                   'raw_pdb': str(raw_path.relative_to(ROOT)), 'scores_json': str(score_path.relative_to(ROOT))}
            audit_rows.append(row)
    writecsv('configuration_consistency.csv', config_rows)
    writecsv('actual_a3m_audit.csv', a3m_rows)
    writecsv('confidence_and_geometry_audit.csv', audit_rows)
    writecsv('recycle_trace.csv', log_rows)
    pkg = Path('/home/g/miniconda3/envs/colabfold_exp0/lib/python3.11/site-packages/colabfold/batch.py')
    lines = pkg.read_text().splitlines()
    selections = [(780, 800), (1355, 1370), (1510, 1528)]
    (OUT / 'installed_colabfold_source_excerpt.txt').write_text(
        '\n'.join([str(pkg), 'sha256=' + sha(pkg)] + [f'{i+1}: {lines[i]}' for lo, hi in selections for i in range(lo-1, hi)]), encoding='utf-8')
    parent = raw_data[20261008, 'BCLW_Bim_A2eT_parent']
    # Save independent chain-level evidence; no model or reward package is loaded.
    writej('input_audit_summary.json', {
        'scope': 'existing artifacts only; no inference, relaxation, scoring, weights loading, or remote MSA',
        'original_files_checked': len(integrity), 'original_file_mismatches': sum(not r['match'] for r in integrity),
        'config_checks': len(config_rows), 'config_mismatches': sum(not r['match'] for r in config_rows),
        'a3m_count': len(a3m_rows), 'a3m_mismatches': sum(not r['exact_query_and_gap_layout'] for r in a3m_rows),
        'structures': len(audit_rows), 'identity_pass': sum(r['identity_raw_and_relaxed'] for r in audit_rows),
        'frozen_gate_pass': sum(r['frozen_gate_pass'] for r in audit_rows),
        'prior_table_disagreements': sum(not r['prior_gate_table_agreement'] for r in audit_rows),
        'max_seq_explanation': 'Installed batch.py lines 1363-1368 clamp single_sequence multimer to effective max_seq=1,max_extra_seq=1; CLI requested 256/512 is recorded correctly.',
        'host_url_explanation': 'Default URL field is configuration metadata, not evidence of a request. Single-sequence branch generates local query rows; no remote MSA call is indicated in logs.',
        'runtime_environment_causality': 'No contemporaneous same-environment positive control; acceleration/backend causes remain untested.',
        'colabfold_batch_path': str(pkg), 'colabfold_batch_sha256': sha(pkg),
        'frozen_direction_criteria_sha256': sha(ROOT / 'output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json')})
    print(json.dumps(readj(OUT / 'input_audit_summary.json'), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
