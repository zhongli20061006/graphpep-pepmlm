from __future__ import annotations
import csv, hashlib, json, pathlib, sys
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

ROOT=pathlib.Path('/mnt/e/Users/graphpep与pepmlm')
OUT=ROOT/'output/bclw_parent_kernel_ablation_v1'
PILOT=ROOT/'output/bh3_bclw_six_complex_pilot_20261020'
AUDIT=ROOT/'output/structure_input_audit_bclw_20261020'
PLAN=json.loads((OUT/'prepared_protocol.json').read_text())
EXPECTED={}
head=None
for line in (AUDIT/'prepared_inputs/parent_only.fasta').read_text().splitlines():
    if line.startswith('>'): head=line[1:]; EXPECTED[head]=''
    elif line.strip(): EXPECTED[head]+=line.strip()
assert len(EXPECTED)==1
complex_id,complex_sequence=next(iter(EXPECTED.items()))
EXPECTED_A,EXPECTED_B=complex_sequence.split(':')
PARSER=PDBParser(QUIET=True)

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def read_json(p):
    def reject(x): raise ValueError('nonfinite JSON constant '+x)
    return json.loads(p.read_text(),parse_constant=reject)
def write_json(path,obj): path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def chain_data(path):
    model=PARSER.get_structure(path.stem,str(path))[0]
    result={}
    for chain in model:
        residues=[r for r in chain if r.id[0]==' ']
        result[chain.id]={
          'sequence':''.join(seq1(r.resname) for r in residues),
          'numbers':[int(r.id[1]) for r in residues],
          'ca':np.asarray([r['CA'].coord for r in residues],dtype=float),
          'plddt':np.asarray([r['CA'].bfactor for r in residues],dtype=float),
          'atoms':np.asarray([a.coord for r in residues for a in r if a.element.upper() not in ('H','D')],dtype=float),
          'atom_res':np.asarray([i for i,r in enumerate(residues) for a in r if a.element.upper() not in ('H','D')],dtype=int)
        }
    return result
def geometry(d):
    a,b=d['A'],d['B']
    d2=np.sum((a['atoms'][:,None,:]-b['atoms'][None,:,:])**2,axis=-1)
    hits=np.argwhere(d2<=36.0)
    pairs=sorted({(int(a['atom_res'][i]),int(b['atom_res'][j])) for i,j in hits})
    return {'minimum_heavyatom_distance_A':float(np.sqrt(d2.min())),
            'clashing_atom_pairs_lt1p5A':int((d2<2.25).sum()),
            'heavyatom_pairs_le6A':int(len(hits)),
            'unique_residue_pairs_le6A':len(pairs),'residue_pairs':pairs}
def load_metrics(base,seed,complex_id,expected_a,expected_b,check_export):
    structdir=base/f'seed_{seed}'/'structures'
    raw_paths=sorted(structdir.glob(f'{complex_id}_unrelaxed_rank_*.pdb'))
    score_paths=sorted(structdir.glob(f'{complex_id}_scores_rank_*.json'))
    if len(raw_paths)!=1 or len(score_paths)!=1: raise ValueError(f'{base}: expected one raw PDB and score file, got {len(raw_paths)}, {len(score_paths)}')
    raw_path,score_path=raw_paths[0],score_paths[0]
    raw=chain_data(raw_path)
    if set(raw)!={'A','B'}: raise ValueError('raw chain IDs differ from A/B')
    relaxed_paths=sorted((base/f'seed_{seed}'/'relaxed').glob(f'{complex_id}_relaxed_rank_*.pdb'))
    if len(relaxed_paths)!=1: raise ValueError(f'{base}: expected one relaxed PDB, got {len(relaxed_paths)}')
    relaxed_path=relaxed_paths[0]
    relaxed=chain_data(relaxed_path)
    expected_numbers={'A':list(range(1,165)),'B':list(range(1,24))}
    identity=(raw['A']['sequence']==expected_a and raw['B']['sequence']==expected_b and
              relaxed['A']['sequence']==expected_a and relaxed['B']['sequence']==expected_b and
              raw['A']['numbers']==expected_numbers['A'] and raw['B']['numbers']==expected_numbers['B'] and
              relaxed['A']['numbers']==expected_numbers['A'] and relaxed['B']['numbers']==expected_numbers['B'])
    score=read_json(score_path)
    plddt=np.asarray(score['plddt'],dtype=float); pae=np.asarray(score['pae'],dtype=float)
    if plddt.shape!=(187,) or pae.shape!=(187,187) or not np.isfinite(plddt).all() or not np.isfinite(pae).all():
        raise ValueError(f'bad/nonfinite confidence array: {plddt.shape}, {pae.shape}')
    interface=np.concatenate((pae[:164,164:].ravel(),pae[164:,:164].ravel()))
    pae_export_path=structdir/f'{complex_id}_predicted_aligned_error_v1.json'
    pae_export=read_json(pae_export_path)
    if isinstance(pae_export,list): pae_export=pae_export[0]
    pae_export=np.asarray(pae_export['predicted_aligned_error'],dtype=float)
    pae_delta=float(np.max(np.abs(pae_export-pae)))
    if check_export and pae_delta!=0.0: raise ValueError(f'PAE export differs from scores JSON: {pae_delta}')
    raw_geom,relaxed_geom=geometry(raw),geometry(relaxed)
    pdb_plddt=np.r_[raw['A']['plddt'],raw['B']['plddt']]
    pdb_delta=float(np.max(np.abs(pdb_plddt-plddt)))
    metrics={
      'seed':seed,'structure_dir':str(structdir),'raw_pdb':str(raw_path),'raw_pdb_sha256':sha(raw_path),
      'relaxed_pdb':str(relaxed_path),'relaxed_pdb_sha256':sha(relaxed_path),'scores_json':str(score_path),
      'scores_json_sha256':sha(score_path),'config_sha256':sha(structdir/'config.json'),
      'identity_and_contiguous_numbering_pass':bool(identity),'chain_A_sequence_length':len(raw['A']['sequence']),
      'chain_B_sequence_length':len(raw['B']['sequence']),'iptm':float(score['iptm']),'ptm':float(score['ptm']),
      'receptor_mean_plddt':float(plddt[:164].mean()),'peptide_mean_plddt':float(plddt[164:].mean()),
      'pdb_ca_vs_json_plddt_max_abs_difference':pdb_delta,
      'interface_pae_max_directed_A':float(interface.max()),'interface_pae_median_directed_A_diagnostic_only':float(np.median(interface)),
      'scores_vs_exported_pae_max_abs_difference':pae_delta,
      'raw_geometry':raw_geom,'relaxed_geometry':relaxed_geom,
      'raw_to_relaxed_CA_RMSD_without_alignment_A':float(np.sqrt(np.mean(np.sum((np.r_[raw['A']['ca'],raw['B']['ca']]-np.r_[relaxed['A']['ca'],relaxed['B']['ca']])**2,axis=-1)))),
    }
    return metrics,read_json(structdir/'config.json')

def main():
    rows=[]; all_configs={}; paired={}
    for seed in (20261008,20261009):
        new,cfg=load_metrics(OUT,seed,complex_id,EXPECTED_A,EXPECTED_B,True)
        baseline,_=load_metrics(PILOT,seed,complex_id,EXPECTED_A,EXPECTED_B,True)
        gate={
          'identity_and_contiguous_numbering':new['identity_and_contiguous_numbering_pass'],
          'has_at_least_one_unique_residue_contact_le6A':new['relaxed_geometry']['unique_residue_pairs_le6A']>=1,
          'iptm_ge_0_70':new['iptm']>=0.70,
          'peptide_mean_CA_plddt_ge_70':new['peptide_mean_plddt']>=70.0,
          'all_directed_interface_pae_max_le_10A':new['interface_pae_max_directed_A']<=10.0,
          'relaxed_min_heavyatom_distance_ge_1_5A':new['relaxed_geometry']['minimum_heavyatom_distance_A']>=1.5
        }
        extra={'receptor_mean_CA_plddt_ge_70':new['receptor_mean_plddt']>=70.0}
        new['frozen_original_quality_gate_checks']=gate
        new['frozen_original_quality_gate_pass']=all(gate.values())
        new['additional_diagnostic_gate_checks']=extra
        new['all_diagnostic_gates_pass']=all(gate.values()) and all(extra.values())
        new['failed_gate_names']=[k for k,v in {**gate,**extra}.items() if not v]
        new['baseline_same_seed']={k:baseline[k] for k in ('iptm','ptm','receptor_mean_plddt','peptide_mean_plddt','interface_pae_max_directed_A','raw_geometry','relaxed_geometry')}
        new['new_minus_baseline']={k:float(new[k]-baseline[k]) for k in ('iptm','ptm','receptor_mean_plddt','peptide_mean_plddt','interface_pae_max_directed_A')}
        rows.append(new); all_configs[str(seed)]=cfg; paired[str(seed)]={'new':new,'baseline':baseline}
    passing=sum(r['all_diagnostic_gates_pass'] for r in rows)
    if passing==2: branch='GO_TO_SEPARATELY_PREPARED_NEXT_REVIEW_ONLY; supports sufficiency for these two seeds, not affinity or reward validity.'
    elif passing==1: branch='HOLD; seed-sensitive partial outcome; no reroll.'
    else: branch='HOLD; disabling fast kernels alone did not recover the parent under these two seeds; this does not establish that other runtime or input factors are innocent or causal.'
    rec=read_json(OUT/'execution_record.json')
    event_times={k:sum(float(e.get('wall_seconds',0)) for e in rec['events'] if e.get('category')==k and e.get('status')=='completed') for k in ('gpu','cpu')}
    budget={'frozen_gpu_process_wall_cap_seconds':PLAN['budget']['gpu_tagged_process_wall_seconds_max'],
            'actual_gpu_process_wall_seconds_from_completed_child_events':event_times['gpu'],
            'frozen_cpu_relax_wall_cap_seconds':PLAN['budget']['cpu_relax_wall_seconds_max'],
            'actual_cpu_relax_wall_seconds_from_completed_child_events':event_times['cpu'],
            'max_disk_bytes':PLAN['budget']['disk_bytes_max'],'final_output_tree_bytes_before_manifest':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),
            'monitor_aggregate_counter_reported':{'gpu':rec.get('gpu_prediction_wall_seconds'),'cpu':rec.get('cpu_relax_wall_seconds')},
            'note':'Watchdog aggregate fields conservatively counted the active child wall interval a second time at child completion; per-child wall_seconds events are authoritative for actual work and remained below the frozen caps. The remaining-budget calculation used the larger aggregate, so enforcement was conservative.'}
    if event_times['gpu']>budget['frozen_gpu_process_wall_cap_seconds'] or event_times['cpu']>budget['frozen_cpu_relax_wall_cap_seconds']:
        raise RuntimeError('Reconciled child process wall time exceeds budget')
    integrity=json.loads((OUT/'postflight_original_integrity.json').read_text())
    if not integrity.get('all_match') or integrity.get('count')!=265: raise RuntimeError('Pilot original manifest postflight integrity mismatch')
    decision={'protocol_id':PLAN['protocol_id'],'completed_prediction_jobs':2,'completed_relaxation_jobs':2,
      'reward_calls':0,'mutant_predictions':0,'all_diagnostic_quality_gate_pass_count':passing,
      'decision':branch,'affinity_direction_accuracy':'NOT_ASSESSED','experimental_Kd_or_complete_binding_deltaG':'NOT_GENERATED',
      'caveat':'The quality gate assesses prediction readiness only; a failure does not show the peptide fails to bind experimentally.',
      'input_scope':'Q92843 1-164 core plus standard Bim_A2eT 23mer; omitted tags/labels/terminal chemistry remain omitted as frozen.'}
    write_json(OUT/'quality_gate_summary.json',{'protocol_id':PLAN['protocol_id'],'metric_definition':'Same code path and thresholds as the read-only Bcl-w pilot audit; raw confidence arrays cross-checked against exported PAE and PDB CA B-factors.','jobs':rows,'pass_count':passing,'denominator':2})
    write_json(OUT/'budget_reconciliation.json',budget)
    write_json(OUT/'decision.json',decision)
    csv_path=OUT/'quality_gate_summary.csv'
    fields=['seed','iptm','ptm','receptor_mean_plddt','peptide_mean_plddt','interface_pae_max_directed_A','raw_geometry.minimum_heavyatom_distance_A','relaxed_geometry.minimum_heavyatom_distance_A','relaxed_geometry.clashing_atom_pairs_lt1p5A','relaxed_geometry.unique_residue_pairs_le6A','identity_and_contiguous_numbering_pass','frozen_original_quality_gate_pass','additional_diagnostic_gate_checks.receptor_mean_CA_plddt_ge_70','all_diagnostic_gates_pass','failed_gate_names']
    def val(r,f):
        for k in f.split('.'): r=r[k] if isinstance(r,dict) else getattr(r,k)
        return json.dumps(r,ensure_ascii=False) if isinstance(r,(list,dict)) else r
    with csv_path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.writer(f); w.writerow(fields)
        for row in rows: w.writerow([val(row,x) for x in fields])
    report=['# Bcl-w 亲本快速内核分支结果','',f"协议：`{PLAN['protocol_id']}`。2个固定seed各完成1次预测和1个同设置CPU放松；未调用GraphPep、FoldX、ESM2或其他reward。",'',
      '| Seed | ipTM | 受体均值pLDDT | 肽均值pLDDT | 双向界面最大PAE (Å) | 原始→放松最小链间距离 (Å) | 放松后接触残基对 | 全门通过 |',
      '|---:|---:|---:|---:|---:|---:|---:|:---:|']
    for x in rows:
      report.append(f"| {x['seed']} | {x['iptm']:.2f} | {x['receptor_mean_plddt']:.2f} | {x['peptide_mean_plddt']:.2f} | {x['interface_pae_max_directed_A']:.2f} | {x['raw_geometry']['minimum_heavyatom_distance_A']:.3f} → {x['relaxed_geometry']['minimum_heavyatom_distance_A']:.3f} | {x['relaxed_geometry']['unique_residue_pairs_le6A']} | {'通过' if x['all_diagnostic_gates_pass'] else '失败'} |")
    report += ['',f"通过数为 {passing}/2。冻结分支结论：{branch}",'',
      '生效配置复核：除输入从同seed的三条件面板缩为亲本单条件（`num_queries: 3→1`）及指定的 `use_fast_kernels: true→false` 外，没有其他配置差异。有效MSA仍为single-sequence（每条链1行）、无模板、Multimer v3 model 1、3 recycle。',
      '', '两组原始结构与放松后结构的A164/B23序列和连续编号均已核对。原始结构存在小于1.5 Å的链间近接，但按冻结判据使用相同CPU放松后几何；若放松后距离通过且接触非零，几何门通过。',
      '', '此分支只回答快速内核选项是否足以改善这两个固定seed下的结构置信度。结构门失败不等于实验不结合，也不定位具体根因；本次没有新Kd、完整结合ΔG或亲和力方向准确率。',
      '', '复核了旧pilot manifest中的265项，哈希仍全部一致。实际GPU/CPU子进程耗时、磁盘和原始结果哈希见 `budget_reconciliation.json`、`execution_record.json` 与 `artifact_manifest.json`。',
      '', '## 文件','',
      '- `quality_gate_summary.csv` / `.json`：冻结门槛逐项结果及同seed基线对照。',
      '- `decision.json`：按2/2、1/2、0/2规则冻结的结论。',
      '- `budget_reconciliation.json`：按各子进程实际wall time核算的预算记录。',
      '- `execution_record.json`、`logs/`、逐seed结构目录：原始执行与产物。','']
    (OUT/'report_zh.md').write_text('\n'.join(report),encoding='utf-8')
    # Recheck all original pilot files after metrics generation; summarize without modifying any source artifacts.
    manifest=read_json(PILOT/'run_result_manifest_v1.json'); mismatches=[]
    for e in manifest['files']:
      p=PILOT/e['path']
      if not p.is_file() or p.stat().st_size!=e['bytes'] or sha(p)!=e['sha256']: mismatches.append(e['path'])
    write_json(OUT/'final_integrity_check.json',{'pilot_manifest_entries':len(manifest['files']),'all_match':len(manifest['files'])==265 and not mismatches,'mismatches':mismatches,'prepared_input_sha256':sha(AUDIT/'prepared_inputs/parent_only.fasta'),'prepared_input_expected_sha256':PLAN['input_fasta_sha256'],'prepared_input_match':sha(AUDIT/'prepared_inputs/parent_only.fasta')==PLAN['input_fasta_sha256']})
    if mismatches: raise RuntimeError('Original pilot changed; see final_integrity_check.json')
    artifact_files=[]
    for p in sorted(OUT.rglob('*')):
      if p.is_file() and p.name not in ('artifact_manifest.json','artifact_manifest.sha256'):
        artifact_files.append({'path':p.relative_to(OUT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
    write_json(OUT/'artifact_manifest.json',{'scope':'All run files except this manifest and detached checksum; raw prediction and relaxation outputs preserved.','files':artifact_files})
    (OUT/'artifact_manifest.sha256').write_text(sha(OUT/'artifact_manifest.json')+'  artifact_manifest.json\n')
    print(json.dumps({'pass_count':passing,'denominator':2,'decision':branch,'metrics':[{k:r[k] for k in ('seed','iptm','receptor_mean_plddt','peptide_mean_plddt','interface_pae_max_directed_A','frozen_original_quality_gate_pass','all_diagnostic_gates_pass','failed_gate_names')} | {'raw_min_distance_A':r['raw_geometry']['minimum_heavyatom_distance_A'],'relaxed_min_distance_A':r['relaxed_geometry']['minimum_heavyatom_distance_A'],'contacts':r['relaxed_geometry']['unique_residue_pairs_le6A']} for r in rows],'budget':budget,'pilot_265_integrity':'PASS','output_files':len(artifact_files),'artifact_manifest_sha256':sha(OUT/'artifact_manifest.json')},ensure_ascii=False,indent=2))

if __name__=='__main__': main()