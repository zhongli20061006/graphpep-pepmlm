"""Paired energy/geometry tables and score-agnostic technical decisions."""
from __future__ import annotations
from collections import defaultdict
import math
from pathlib import Path
import statistics
from pipeline.foldx_matched_diagnostic import require, save_json, sha, csv_text, read_csv, summarize_backbone, summarize_seed
from pipeline import backbone_geometry as geo

COUNT_TERMS=frozenset(('Number of Residues','Interface Residues','Interface Residues Clashing','Interface Residues VdW Clashing','Interface Residues BB Clashing'))


def term_units(name: str) -> str:
    return 'count' if name in COUNT_TERMS else 'kcal/mol'


def write_csv(path: Path, rows: list[dict], fields=None):
    if fields is None:
        fields=list(dict.fromkeys(k for row in rows for k in row))
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(csv_text(rows,fields),encoding='utf-8',newline='\n')


def paired_contrasts(rows: list[dict]) -> list[dict]:
    groups=defaultdict(dict)
    for row in rows:
        require(row['condition'] not in groups[row['pair_uid']],'Duplicate corresponding condition')
        groups[row['pair_uid']][row['condition']]=row
    result=[]
    for uid,conditions in sorted(groups.items()):
        first=next(iter(conditions.values()))
        row={k:first[k] for k in ('pair_uid','peptide_id','seed','model_number','repeat')}
        complete=set(conditions)=={'fixed','flexible'} and all(v['status']=='ok' for v in conditions.values())
        row.update(status='ok' if complete else 'missing',failure_reason='' if complete else 'Missing/failed corresponding fixed or flexible pair',fixed_gain_kcal_mol=conditions.get('fixed',{}).get('gain_kcal_mol',''),flexible_gain_kcal_mol=conditions.get('flexible',{}).get('gain_kcal_mol',''),flexible_minus_fixed_gain_kcal_mol=conditions['flexible']['gain_kcal_mol']-conditions['fixed']['gain_kcal_mol'] if complete else '')
        result.append(row)
    return result


def contrast_summary(rows: list[dict]) -> list[dict]:
    groups=defaultdict(list)
    for row in rows: groups[(row['peptide_id'],row['seed'],row['model_number'])].append(row)
    results=[]
    for key,items in sorted(groups.items()):
        values=[x['flexible_minus_fixed_gain_kcal_mol'] for x in items if x['status']=='ok']
        complete=len(items)==5 and len(values)==5 and {int(x['repeat']) for x in items}==set(range(1,6))
        results.append(dict(zip(('peptide_id','seed','model_number'),key),status='complete' if complete else 'incomplete',n_complete_contrasts=len(values),median_flexible_minus_fixed_gain_kcal_mol=statistics.median(values) if complete else '',minimum_contrast_kcal_mol=min(values) if values else '',maximum_contrast_kcal_mol=max(values) if values else ''))
    return results


def technical_gate(rows: list[dict], expected: int) -> dict:
    failures=[]
    if len(rows)!=expected or len({x['job_id'] for x in rows})!=expected: failures.append('Incomplete/duplicate expected structure jobs')
    if any(x['status']!='ok' for x in rows): failures.append('Missing/failed jobs')
    valid=[x for x in rows if x['status']=='ok']
    if any(not math.isfinite(float(x['energy_kcal_mol'])) for x in valid): failures.append('Nonfinite actual Interaction Energy')
    if any(float(x.get('frozen_max_A',math.inf))>1e-6 for x in valid): failures.append('Frozen atom displacement exceeds tolerance')
    if not any(x['condition']=='flexible' and float(x.get('B_backbone_max_A',0))>.01 for x in valid): failures.append('No measurable flexible peptide backbone movement >0.01 Angstrom')
    return dict(technical_pass=not failures,expected_structures=expected,observed_structures=len(rows),valid_structures=len(valid),converged_structures=sum(x.get('converged',False) for x in valid),nonconverged_structures=sum(not x.get('converged',False) for x in valid),failure_reasons=failures,gate_uses_score_direction=False)


def assemble_pairs(job_manifest: list[dict], results: dict) -> tuple[list[dict],list[dict],list[dict]]:
    structures=[]; terms=[]; grouped=defaultdict(dict)
    for job in job_manifest:
        result=results.get(job['job_id'],{})
        row={k:job[k] for k in ('job_id','source_uid','pair_uid','peptide_id','seed','model_number','repeat','role','condition','experimental_direction')}
        for key in ('seed','model_number','repeat'): row[key]=int(row[key])
        row.update(status=result.get('status','missing'),failure_reason=result.get('failure_reason','Job not run'),provenance_valid=result.get('provenance_valid',False),energy_kcal_mol='',units='kcal/mol')
        if row['status']=='ok':
            row['energy_kcal_mol']=result['raw']['energy_kcal_mol']; row['failure_reason']=''
            relaxation=result['relaxation']
            row.update({k:v for k,v in relaxation.items() if k in ('frozen_max_A','converged','iterations','iteration_cap','iteration_cap_reached','hydration_sha256','start_positions_sha256') or k.startswith(('A_backbone_','B_backbone_','B_all_','B8_12_','alpha_like_','hbond_'))})
            for time in ('before','after'):
                row.update({time+'_'+k:v for k,v in relaxation[time].items()})
            for table,values in result['raw']['tables'].items():
                if not table.startswith('Interaction_'): continue
                for name,value in values.items(): terms.append({k:row[k] for k in ('job_id','source_uid','pair_uid','peptide_id','seed','model_number','repeat','role','condition')}|dict(table_path=f"jobs/{job['job_id']}/{table}",term=name,value=value,units=term_units(name)))
        structures.append(row)
        group=grouped[(job['pair_uid'],job['condition'])]
        require(job['role'] not in group,'Duplicate source pair role')
        group[job['role']]=row
    pairs=[]
    for (uid,condition),roles in sorted(grouped.items()):
        first=next(iter(roles.values()))
        row={k:first[k] for k in ('pair_uid','peptide_id','seed','model_number','repeat','condition','experimental_direction')}
        ok=set(roles)=={'wt_reference','mutant'} and all(x['status']=='ok' for x in roles.values())
        row.update(status='ok' if ok else 'missing',failure_reason='' if ok else '; '.join(role+': '+str(roles.get(role,{}).get('failure_reason','role missing')) for role in ('wt_reference','mutant') if roles.get(role,{}).get('status')!='ok'),wt_job_id=roles.get('wt_reference',{}).get('job_id',''),mutant_job_id=roles.get('mutant',{}).get('job_id',''),wt_energy_kcal_mol=roles.get('wt_reference',{}).get('energy_kcal_mol',''),mutant_energy_kcal_mol=roles.get('mutant',{}).get('energy_kcal_mol',''),gain_kcal_mol=roles['wt_reference']['energy_kcal_mol']-roles['mutant']['energy_kcal_mol'] if ok else '',units='kcal/mol')
        pairs.append(row)
    return structures,pairs,terms


def summaries(pairs: list[dict], contrasts: list[dict]):
    groups=defaultdict(list)
    for row in pairs: groups[(row['peptide_id'],row['condition'],row['seed'],row['model_number'],row['experimental_direction'])].append(row)
    backbones=[]
    cgroups={(x['peptide_id'],x['seed'],x['model_number']):x for x in contrast_summary(contrasts)}
    for key,items in sorted(groups.items()):
        summary=summarize_backbone([x['gain_kcal_mol'] for x in items if x['status']=='ok'])
        summary.update(dict(zip(('peptide_id','condition','seed','model_number','experimental_direction'),key)))
        summary['median_flexible_minus_fixed_gain_kcal_mol']=cgroups[(key[0],key[2],key[3])]['median_flexible_minus_fixed_gain_kcal_mol']
        backbones.append(summary)
    seeds=[]; sg=defaultdict(list)
    for row in backbones: sg[(row['peptide_id'],row['condition'],row['seed'],row['experimental_direction'])].append(row)
    for key,items in sorted(sg.items()):
        summary=summarize_seed([x['median_gain_kcal_mol'] for x in items if x['complete']],key[3])
        summary.update(dict(zip(('peptide_id','condition','seed','experimental_direction'),key)))
        cs=[cgroups[(key[0],key[2],model)]['median_flexible_minus_fixed_gain_kcal_mol'] for model in (1,2,3)]
        summary['median_flexible_minus_fixed_gain_kcal_mol']=statistics.median(cs) if all(isinstance(x,(int,float)) for x in cs) else ''
        seeds.append(summary)
    return backbones,seeds


def publish_hashes(output: Path):
    analysis=output/'analysis'
    owned=('structure_results.csv','pair_results.csv','contrasts.csv','backbone_summary.csv','seed_summary.csv','energy_terms.csv','torsions.csv','hbond_proxies.csv','residue_displacements.csv','pilot_decision.json','decision.json','report_zh.md','baseline_geometry.csv','baseline_torsions.csv','baseline_hbond_proxies.csv','baseline_energy_term_contrasts.csv','baseline.json')
    save_json(analysis/'published_artifact_hashes.json',{(analysis/name).relative_to(output).as_posix():sha(analysis/name) for name in owned if (analysis/name).is_file()})


def analyze(output: Path, jobs: list[dict], results: dict, phase='full') -> dict:
    analysis=output/'analysis'; analysis.mkdir(parents=True,exist_ok=True)
    structures,pairs,terms=assemble_pairs(jobs,results); contrasts=paired_contrasts(pairs); backbones,seeds=summaries(pairs,contrasts)
    for name,rows in [('structure_results',structures),('pair_results',pairs),('contrasts',contrasts),('backbone_summary',backbones),('seed_summary',seeds),('energy_terms',terms)]: write_csv(analysis/(name+'.csv'),rows,fields=None if rows else ['job_id','term','value','units'])
    detail={'torsions':[],'hbond_proxies':[],'residue_displacements':[]}
    for job in jobs:
        result=results.get(job['job_id'],{})
        if result.get('status')!='ok': continue
        meta={k:job[k] for k in ('job_id','source_uid','pair_uid','peptide_id','seed','model_number','repeat','role','condition')}
        relaxation=result['relaxation']
        for moment in ('before','after'):
            detail['torsions'].extend(meta|dict(moment=moment)|row for row in relaxation['torsions_'+moment])
            detail['hbond_proxies'].extend(meta|dict(moment=moment)|row for row in relaxation['hbond_proxies_'+moment])
        detail['residue_displacements'].extend(meta|row for row in relaxation['residue_displacements'])
    for name,rows in detail.items(): write_csv(analysis/(name+'.csv'),rows,fields=None if rows else ['job_id','source_uid','moment','residue'])
    pilot=[x for x in structures if int(x['repeat'])==1]
    pilot_decision=technical_gate(pilot,48)
    pilot_decision.update(phase='pilot',expected_pairs=24,complete_pairs=sum(x['status']=='ok' and int(x['repeat'])==1 for x in pairs),expand_to_full=pilot_decision['technical_pass'])
    save_json(analysis/'pilot_decision.json',pilot_decision)
    decision=technical_gate(structures,240)
    decision.update(phase='full',complete=len(structures)==240 and all(x['status']=='ok' for x in structures) and len(pairs)==120 and all(x['status']=='ok' for x in pairs),expected_pairs=120,complete_pairs=sum(x['status']=='ok' for x in pairs),expected_contrasts=60,complete_contrasts=sum(x['status']=='ok' for x in contrasts),diagnostic_directions=[{k:x[k] for k in ('peptide_id','condition','seed','diagnostic_direction')} for x in seeds],biological_replicates=False)
    save_json(analysis/'decision.json',decision)
    lines=['# 固定与可动肽骨架开发诊断','',f"技术试点：{pilot_decision['valid_structures']}/48 个结构有效；技术门槛{'通过' if pilot_decision['technical_pass'] else '未通过'}。",f"完整批次：{decision['valid_structures']}/240 个结构、{decision['complete_pairs']}/120 个对应配对、{decision['complete_contrasts']}/60 个条件差值有效；完整状态：{decision['complete']}。",f"收敛报告：{decision['converged_structures']} 个达到可动原子力 RMS 阈值，{decision['nonconverged_structures']} 个未达到；迭代上限不等于收敛。",'', '结合分数仅为实际 FoldX AnalyseComplex 的 A/B Interaction Energy。增益 = 对应 WT reference − mutant；条件差值先在同一次重复内计算 flexible − fixed，再取五次重复中位数和三个模型中位数。技术门槛不使用生物学方向。','', '|突变|条件|种子|方向判定|中位增益 kcal/mol|','|---|---|---|---|---|']
    lines += [f"|{x['peptide_id']}|{x['condition']}|{x['seed']}|{x['diagnostic_direction']}|{x['median_gain_kcal_mol']}|" for x in seeds]
    lines += ['', '这是筛选后的开发敏感性诊断，不是留出验证。局部隐式溶剂最小化不包含热采样、未结合态或新 Kd；弱坐标约束没有指定螺旋。alpha-like 扭转角区间及 O_i--N_i+4 距离仅为描述性几何代理，不是 DSSP 或氢键占有率。重复和 AFM 模型不是生物学重复。分项能量不要求和为结合能；PDB 舍入与力场误差仍存在。','']
    (analysis/'report_zh.md').write_text('\n'.join(lines),encoding='utf-8')
    publish_hashes(output)
    return pilot_decision if phase=='pilot' else decision


def baseline(output: Path, root: Path, protocol: dict) -> dict:
    analysis=output/'analysis'; analysis.mkdir(parents=True,exist_ok=True)
    structures=[]; torsions=[]; hbond=[]; cached={}
    from pipeline.foldx_matched_diagnostic import validate_pdb
    sequences={'P53_WT':protocol['scientific_protocol']['parent_sequence'],**{m['peptide_id']:m['sequence'] for m in protocol['scientific_protocol']['mutations']}}
    for seed in protocol['scientific_protocol']['seeds']:
        for model in protocol['scientific_protocol']['models']:
            for pid,seq in sequences.items():
                name=f'{pid}_seed{seed}_model{model}'
                path=Path(protocol['independent_afm_root'])/'jobs'/name/'input_Repair.pdb'
                text=(root/path).read_text(); validate_pdb(text,seq); coords=geo.coordinates(text); cached[pid,seed,model]=coords
                meta=dict(peptide_id=pid,seed=seed,model_number=model,source_path=path.as_posix(),source_sha256=sha(root/path))
                torsions.extend(meta|x for x in geo.torsions(coords,seq)); hbond.extend(meta|x for x in geo.hbond_proxies(coords,seq))
                aligned=geo.align_on_receptor(cached['P53_WT',seed,model],coords)
                reference=cached['P53_WT',seed,model]
                # Mutations change sidechain atom identities; compare corresponding backbone only.
                row=meta|geo.geometry(coords,seq)
                for chain,bb,start,end,prefix in [('A',True,None,None,'A_backbone'),('B',True,None,None,'B_backbone'),('B',True,8,12,'B8_12_backbone')]:
                    row.update({prefix+'_'+k:v for k,v in geo.displacement(geo.select(reference,chain,bb,start,end),geo.select(aligned,chain,bb,start,end)).items()})
                row['alignment']='receptor A CA Kabsch'; structures.append(row)
    old=read_csv(root/protocol['source_energy_terms_csv']); groups=defaultdict(dict)
    for row in old:
        if row['peptide_id'] not in sequences or row['peptide_id']=='P53_WT' or not Path(row['table_path']).name.startswith('Interaction_'): continue
        key=(row['job_id'],row['peptide_id'],int(row['seed']),int(row['model_number']),int(row['repeat']),row['term'])
        require(row['role'] not in groups[key],'Duplicate original Interaction term')
        groups[key][row['role']]=float(row['value'])
    terms=[]
    for key,roles in sorted(groups.items()):
        require(set(roles)=={'wt_reference','mutant'},'Missing old matched energy term role')
        terms.append(dict(zip(('job_id','peptide_id','seed','model_number','repeat','term'),key),wt_reference=roles['wt_reference'],mutant=roles['mutant'],wt_minus_mutant=roles['wt_reference']-roles['mutant'],units=term_units(key[-1])))
    require(len({(x['job_id'],x['repeat']) for x in terms})==60,'Expected 60 original matched term pairs')
    for name,rows in [('baseline_geometry',structures),('baseline_torsions',torsions),('baseline_hbond_proxies',hbond),('baseline_energy_term_contrasts',terms)]: write_csv(analysis/(name+'.csv'),rows)
    report=dict(independent_structures=len(structures),original_matched_pairs=60,alignment='A CA Kabsch only for independently folded structures',summary_tables_excluded=True,energy_term_sum_asserted=False)
    save_json(analysis/'baseline.json',report); publish_hashes(output)
    return report
