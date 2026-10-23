"""Map deposited sequences/observations to the assay core without generating poses."""
from __future__ import annotations

import json
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from audit_inputs import OUT, ROOT, SOURCE, writecsv, writej


def main() -> None:
    inventory = json.loads((OUT / 'structure_database_inventory.json').read_text())
    receptor = ''.join((SOURCE / 'raw/Q92843_BCLW.fasta').read_text().splitlines()[1:])[:164]
    parent = 'RPEIWITQELRRIGDEFNAYYAR'
    rows, difference_rows = [], []
    summaries = {}
    for entry in inventory['entries']:
        for entity in entry['entities']:
            seq = entity['sequence']
            # Assignment uses explicit entry/entity identity plus sequence, not ambiguous molecule labels.
            is_bclw = ((entry['pdb_id'] != '3FDL' and entity['entity_id'] == '1' and entry['pdb_id'] != '4K5B') or
                       (entry['pdb_id'] == '4K5B' and entity['entity_id'] == '2'))
            is_bim = entry['pdb_id'] == '3FDL' and entity['entity_id'] == '2'
            target = receptor if is_bclw else parent if is_bim else None
            row = {'pdb_id': entry['pdb_id'], 'entity': entity['entity_id'], 'author_chains': entity['chains'],
                   'deposited_description': entity['description'], 'deposited_length': len(seq),
                   'deposited_mutation_annotation': entity['mutation'],
                   'assigned_role': 'Bcl-w related receptor' if is_bclw else 'Bim related peptide' if is_bim else 'other partner',
                   'target_mapping': 'Q92843 1-164' if is_bclw else 'Bim_A2eT 23mer' if is_bim else 'not aligned to assay component',
                   'target_positions_covered': '', 'identical_positions': '', 'substitutions_relative_to_current_input': '',
                   'unmapped_target_positions': '', 'extra_deposited_residues': ''}
            if target:
                # Audited entries have contiguous cores without internal indels.
                # Use a unique conserved N-terminal anchor to exclude affinity tags;
                # do not match target Met1 to an unrelated Met/Ser inside the tag.
                if is_bclw:
                    anchor = receptor[:28]
                    start = seq.find(anchor)
                    target_start = 0
                    if start < 0:
                        anchor = receptor[1:28]
                        start = seq.find(anchor)
                        target_start = 1
                    assert start >= 0 and seq.count(anchor) == 1
                else:
                    start, target_start = seq.index('RPEIW'), 0
                count = min(len(target)-target_start, len(seq)-start)
                mapped = {target_start+i: start+i for i in range(count)}
                differences = [f'{target[i]}{i+1}{seq[j]}' for i, j in mapped.items() if target[i] != seq[j]]
                missing = [i+1 for i in range(len(target)) if i not in mapped]
                extra = [j+1 for j in range(len(seq)) if j not in mapped.values()]
                row.update(target_positions_covered=len(mapped), identical_positions=sum(target[i] == seq[j] for i,j in mapped.items()),
                           substitutions_relative_to_current_input=';'.join(differences),
                           unmapped_target_positions=';'.join(map(str, missing)), extra_deposited_residues=';'.join(map(str, extra)))
                difference_rows.extend({'pdb_id': entry['pdb_id'], 'entity': entity['entity_id'],
                                        'target_position': i+1, 'target_residue': target[i], 'entity_seq_position': j+1,
                                        'deposited_residue': seq[j], 'identical': target[i] == seq[j]} for i,j in mapped.items())
            rows.append(row)
        path = OUT / 'sources' / (entry['pdb_id'] + '.cif')
        if path.exists():
            cif = MMCIF2Dict(str(path))
            model_ids = cif['_atom_site.pdbx_PDB_model_num']
            model1 = min(map(int, model_ids))
            seen = {}
            for mod, asym, auth, entity_id, number, atom in zip(model_ids, cif['_atom_site.label_asym_id'],
                    cif['_atom_site.auth_asym_id'], cif['_atom_site.label_entity_id'],
                    cif['_atom_site.label_seq_id'], cif['_atom_site.label_atom_id']):
                if int(mod) == model1 and atom == 'CA' and number not in ('?', '.'):
                    key = f'{asym}/auth:{auth}/entity:{entity_id}'
                    seen.setdefault(key, set()).add(int(number))
            summaries[entry['pdb_id']] = {'model_count': len(set(model_ids)),
                'method': entry['method'], 'observed_ca_count_first_model': {key:len(value) for key,value in seen.items()},
                'observed_label_seq_ids_first_model': {key:sorted(value) for key,value in seen.items()},
                'struct_ref_seq_dif': {key:value for key,value in cif.items() if key.startswith('_struct_ref_seq_dif.')},
                'mapping_note': 'Observed counts are first model only; PDB label_seq_id maps to deposited entity sequence, not automatically UniProt numbering. Sequence correspondence uses unique conserved N-terminal anchors and contiguous core mapping; tags are excluded.'}
    writecsv('deposited_sequence_mapping.csv', rows)
    writecsv('residue_correspondence.csv', difference_rows)
    writej('coordinate_observation_audit.json', summaries)
    # Compare the historical working setup without transferring its runtime validation.
    historical = ROOT / 'output/feedback_control_1ycr_20261006/structures/config.json'
    old = json.loads(historical.read_text())
    current = json.loads((ROOT / 'output/bh3_bclw_six_complex_pilot_20261020/seed_20261008/structures/config.json').read_text())
    keys = ['version','msa_mode','max_seq','max_extra_seq','model_type','model_order','num_models','num_recycles',
            'use_fast_kernels','use_bfloat16','use_fuse','compile_mode','use_templates']
    writej('historical_positive_setup_comparison.json', {'historical_path': str(historical.relative_to(ROOT)),
           'current_path': 'output/bh3_bclw_six_complex_pilot_20261020/seed_20261008/structures/config.json',
           'settings': [{'field':key,'historical':old.get(key),'current':current.get(key)} for key in keys],
           'interpretation':'Different target and MSA/input data; same version/options do not establish a contemporaneous controlled positive or prove the current implementation is numerically correct.'})
    for r in rows:
        print(r['pdb_id'], r['entity'], r['assigned_role'], r['substitutions_relative_to_current_input'], 'missing',r['unmapped_target_positions'])


if __name__ == '__main__':
    main()
