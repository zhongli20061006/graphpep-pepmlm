"""Verify scientific audit cross-links, numeric claims, immutable inputs and prepared-only scope."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

from audit_inputs import OUT, PILOT, ROOT, sha, writej


def main() -> None:
    checks: list[dict] = []
    def check(name: str, ok: bool, detail: object = None) -> None:
        checks.append({'check':name,'pass':bool(ok),'detail':detail})
    def load(name: str) -> object:
        return json.loads((OUT/name).read_text(encoding='utf-8-sig'))
    def rows(name: str) -> list[dict]:
        with (OUT/name).open(encoding='utf-8-sig') as stream:
            return list(csv.DictReader(stream))
    manifest=json.loads((PILOT/'run_result_manifest_v1.json').read_text())
    errors=[]
    for r in manifest['files']:
        path=PILOT/r['path']
        if not path.exists() or path.stat().st_size!=r['bytes'] or sha(path)!=r['sha256']:
            errors.append(r['path'])
    check('original_265_files_still_unchanged',len(manifest['files'])==265 and not errors,errors)
    expected='23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b'
    check('direction_criteria_unchanged',sha(ROOT/'output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json')==expected)
    # Parse every delivered JSON/CSV; verify row widths independently of generation.
    parsed=[]
    for p in OUT.rglob('*.json'):
        if p.name not in ('verification.json','artifact_manifest.json'):
            json.loads(p.read_text(encoding='utf-8-sig'));parsed.append(str(p.relative_to(OUT)))
    for p in OUT.rglob('*.csv'):
        with p.open(encoding='utf-8-sig') as stream:
            data=list(csv.reader(stream))
        check('csv_row_width:'+p.name, bool(data) and all(len(r)==len(data[0]) for r in data))
    check('all_json_parse',True,parsed)
    for p in OUT.glob('*.py'):
        compile(p.read_text(),str(p),'exec')
    check('audit_scripts_compile_without_inference',True)
    a=rows('confidence_and_geometry_audit.csv')
    check('six_rows_keep_failure_denominator',len(a)==6 and all(r['frozen_gate_pass']=='False' for r in a))
    check('all_identity_and_old_table_agreement',all(r['identity_raw_and_relaxed']=='True' and r['prior_gate_table_agreement']=='True' for r in a))
    check('confidence_exports_exact',all(float(r['exported_PAE_vs_scores_PAE_max_abs_difference'])==0 and float(r['pdb_ca_vs_json_max_abs_plddt_difference'])==0 for r in a))
    rec=[float(r['receptor_mean_plddt']) for r in a]
    check('reported_receptor_means',abs(min(rec)-39.91823170731708)<1e-9 and abs(max(rec)-41.32384146341463)<1e-9,
          {'minimum':min(rec),'maximum':max(rec)})
    check('no_receptor_residue_ge70',all(float(r['receptor_fraction_plddt_ge70'])==0 for r in a))
    check('msa_queries_only',len(rows('actual_a3m_audit.csv'))==6 and all(r['exact_query_and_gap_layout']=='True' and r['homolog_sequences']=='0' for r in rows('actual_a3m_audit.csv')))
    check('34_effective_setting_checks',len(rows('configuration_consistency.csv'))==34 and all(r['match']=='True' for r in rows('configuration_consistency.csv')))
    check('all_recycle_steps_retained',len(rows('recycle_trace.csv'))==24 and
          all({int(r['recycle']) for r in rows('recycle_trace.csv') if r['seed']==str(seed) and r['complex_id']==name}=={0,1,2,3}
              for seed in (20261008,20261009) for name in {r['complex_id'] for r in a}))
    inventory=load('structure_database_inventory.json')
    check('eight_entry_bounded_structure_review',len(inventory['entries'])==8 and inventory['total_query_count']==5 and len(rows('structure_source_input_map.csv'))==8)
    mapping=rows('deposited_sequence_mapping.csv')
    check('1ZY3_current_numbering_differences',next(r for r in mapping if r['pdb_id']=='1ZY3' and r['entity']=='1')['substitutions_relative_to_current_input']=='P117V;Q133R')
    check('4K5B_role_resolved_by_sequence',next(r for r in mapping if r['pdb_id']=='4K5B' and r['entity']=='2')['author_chains']=='C,D')
    uni=load('sources/Q92843_current.json')
    variants=[f for f in uni['features'] if f.get('featureId')=='VAR_048418']
    check('Q133R_natural_variant_primary_record',len(variants)==1 and variants[0]['location']['start']['value']==133 and
          variants[0]['alternativeSequence']=={'originalSequence':'Q','alternativeSequences':['R']})
    source_errors=[]
    for row in load('source_retrieval_manifest.json'):
        if row.get('path'):
            path=ROOT/row['path']
            if not path.exists() or path.stat().st_size!=row['bytes'] or sha(path)!=row['sha256']:
                source_errors.append(row['path'])
    check('all_reused_and_new_primary_source_hashes',not source_errors,source_errors)
    plan=load('next_compute_protocol_prepared_v1.json')
    check('prepared_input_hash',sha(Path(plan['input_fasta']))==plan['input_fasta_sha256'])
    original_events=json.loads((PILOT/'run_execution_record_v1.json').read_text())['events']
    for cmd in plan['commands']:
        old=list(next(e for e in original_events if e['label']==f"afm_seed_{cmd['seed']}")['command'])
        expected_argv=old.copy()
        expected_argv[1]=plan['input_fasta']
        expected_argv[2]=str(Path(plan['planned_output_root'])/f"seed_{cmd['seed']}/structures")
        expected_argv[expected_argv.index('--use-fast-kernels')]='--no-use-fast-kernels'
        check('planned_only_kernel_option_changed:'+str(cmd['seed']),cmd['argv']==expected_argv)
    old_config=json.loads((PILOT/'run_config_frozen_v1.json').read_text())
    check('old_quality_gates_kept',plan['quality_gate_unchanged']==old_config['proposed_readouts_pre_score']['quality_gate'])
    check('extra_receptor_gate_tightens_only',plan['additional_diagnostic_gate']['receptor_mean_plddt_min']==70)
    check('no_new_run_directory_exists',not Path(plan['planned_output_root']).exists())
    decision=load('decision.json')
    check('no_new_inference_or_reward_calls_reported',decision['new_prediction_jobs_started']==0 and decision['new_reward_calls']==0 and 'NOT_STARTED' in plan['status'])
    check('bounded_proposal',len(plan['commands'])==2 and plan['budget']['new_prediction_jobs_max']==2 and plan['budget']['no_reward_calls'])
    missing_links=[]
    for p in OUT.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if '://' not in link and not (p.parent/link.split('#')[0]).exists():
                missing_links.append({'file':p.name,'link':link})
    check('markdown_local_links_resolve',not missing_links,missing_links)
    check('hypotheses_and_assay_tables_present',len(rows('root_cause_hypotheses_evidence.csv'))==8 and len(rows('assay_to_input_correspondence.csv'))==10)
    result={'date':'2026-10-20 Asia/Shanghai','verification_scope':'file identity, numeric cross-checks, source mappings and prepared-only protocol; not model or biological validity',
            'check_count':len(checks),'failed_checks':sum(not c['pass'] for c in checks),'checks':checks}
    writej('verification.json',result)
    if result['failed_checks']:
        print(json.dumps([c for c in checks if not c['pass']],ensure_ascii=False,indent=2))
        raise SystemExit(1)
    files=[]
    for path in sorted(OUT.rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.name not in ('artifact_manifest.json','artifact_manifest.sha256'):
            files.append({'path':str(path.relative_to(OUT)),'bytes':path.stat().st_size,'sha256':sha(path)})
    writej('artifact_manifest.json',{'scope':'All non-ephemeral audit files except manifest and detached checksum; original input/source references verified separately.',
                                  'files':files})
    (OUT/'artifact_manifest.sha256').write_text(sha(OUT/'artifact_manifest.json')+'  artifact_manifest.json\n')
    print(json.dumps({'checks_passed':len(checks),'failed_checks':0,'original_files_unchanged':265,
                      'audit_manifest_files':len(files),'artifact_manifest_sha256':sha(OUT/'artifact_manifest.json')},indent=2))


if __name__ == '__main__':
    main()
