"""Frozen fixed-versus-flexible paired backbone diagnostic.

Prepare/baseline/run/analyze/verify are separate auditable commands. Every score
is reparsed from actual AnalyseComplex tables. Full results require 240 valid
structure jobs; failures remain explicit and cannot be promoted on resume.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

if __package__ in (None,''):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pipeline.foldx_matched_diagnostic import require, sha, save_json, load_json, read_csv, csv_text, check_hashes, validate_pdb, raw_energy, command_log, RECEPTOR, SETTINGS
from pipeline import backbone_geometry as geo
from pipeline import backbone_analysis as analysis

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'output/backbone_diagnostic_20261007'
CODE=('pipeline/backbone_diagnostic.py','pipeline/backbone_relax.py','pipeline/backbone_geometry.py','pipeline/backbone_analysis.py','tests/test_backbone_diagnostic.py')
MUTATIONS={'P53_P27A':'ETFSDLWKLLAE','P53_L22A':'ETFSDAWKLLPE'}
PARENT='ETFSDLWKLLPE'
ANALYSE_CONFIG='command=AnalyseComplex\npdb=relaxed.pdb\nanalyseComplexChains=A,B\n'+SETTINGS


def analyse_command(binary: Path):
    return [str(binary),'-f','analyse.cfg']


def validate_source(text: str, sequence: str) -> dict:
    require(sequence in (PARENT,*MUTATIONS.values()),'Unsupported mutation sequence')
    validate_pdb(text,sequence)
    return geo.coordinates(text)


def parse_score(directory: Path, name: str) -> dict:
    result=raw_energy(directory,name)
    for terms in result['tables'].values():
        require(all(k in terms and math.isfinite(terms[k]) for k in ('Interaction Energy','IntraclashesGroup1','IntraclashesGroup2')),'Missing/nonfinite actual intraclashes')
    return result


def implementation_hashes(root=ROOT):
    return {name:sha(root/name) for name in CODE}


def validate_protocol(p: dict):
    require(p['state']=='FROZEN_BEFORE_NEW_BACKBONE_RELAXATION_AND_FOLDX_SCORES','Protocol not frozen before new scores')
    s=p['scientific_protocol']
    expected=dict(receptor_sequence=RECEPTOR,parent_sequence=PARENT,seeds=[20261011,20261012],models=[1,2,3],repeats=[1,2,3,4,5],pilot_repeats=[1],roles=['wt_reference','mutant'],conditions=['fixed','flexible'],chains=['A','B'],free_N_C_termini=True,backbone_atoms=['N','CA','C','O','OXT'],peptide_backbone_tether_k_kcal_mol_A2=.1,peptide_backbone_tether_k_kJ_mol_nm2=41.84,forcefields=['amber14/protein.ff14SB.xml','implicit/obc2.xml'],solute_dielectric=1,solvent_dielectric=78.5,implicit_solvent_kappa_per_nm=1.274,nonbonded_method='NoCutoff',constraints=None,rigidWater=False,removeCMMotion=False,pH=7.4,hydrogen_seed_base=20261007,openmm_platform='CPU',openmm_threads=1,openmm_deterministic_forces=True,minimization_tolerance_kJ_mol_nm=10,max_minimization_iterations=2000,foldx_command='AnalyseComplex',foldx_chains='A,B',foldx_temperature_K=298,foldx_pH=7.4,foldx_ionStrength_M=.15,foldx_water='-IGNORE',foldx_noNterm=False,foldx_noCterm=False,pilot_expected_structures=48,pilot_expected_pairs=24,full_expected_structures=240,full_expected_pairs=120,frozen_atom_tolerance_A=1e-6)
    require(all(s.get(k)==v for k,v in expected.items()),'Frozen scientific protocol/settings mismatch')
    require(s['mutations']==[dict(peptide_id='P53_P27A',sequence=MUTATIONS['P53_P27A'],position=11,experimental_direction='stronger'),dict(peptide_id='P53_L22A',sequence=MUTATIONS['P53_L22A'],position=6,experimental_direction='weaker')],'Frozen mutation identity mismatch')


def verify_environment(protocol: dict):
    environment=protocol['environment']
    check_hashes(environment['external_hashes'],Path('/'))
    import openmm
    require(openmm.__version__==environment['openmm_version'],'Frozen OpenMM version changed')
    require(Path(sys.executable).resolve()==Path(environment['python_executable']).resolve(),'Use the frozen WSL Python executable')


def build_manifests(protocol: dict, root: Path):
    pairs=[x for x in read_csv(root/protocol['source_pair_csv']) if x['peptide_id'] in MUTATIONS]
    expected={(pid,seed,model,repeat) for pid in MUTATIONS for seed in (20261011,20261012) for model in (1,2,3) for repeat in range(1,6)}
    observed={(x['peptide_id'],int(x['seed']),int(x['model_number']),int(x['repeat'])) for x in pairs}
    require(len(pairs)==60 and observed==expected,'Missing/duplicate frozen source pair identities')
    sources=[]; jobs=[]
    oldroot=Path(protocol['source_pair_csv']).parent.parent
    for pair in sorted(pairs,key=lambda x:(x['peptide_id'],int(x['seed']),int(x['model_number']),int(x['repeat']))):
        require(pair['status']=='ok' and pair['provenance_valid']=='True','Old source pair not valid')
        pid=pair['peptide_id']; seed=int(pair['seed']); model=int(pair['model_number']); repeat=int(pair['repeat'])
        oldjob=f'{pid}_seed{seed}_model{model}'
        require(pair['job_id']==oldjob and pair['parent_id']=='P53_WT','Source pair metadata identity mismatch')
        direction='stronger' if pid=='P53_P27A' else 'weaker'
        require(pair['experimental_direction']==direction,'Source experimental direction mismatch')
        uid=f'{oldjob}_repeat{repeat}'
        for role,column,sequence in [('wt_reference','wt_path',PARENT),('mutant','mutant_path',MUTATIONS[pid])]:
            path=oldroot/Path(pair[column])
            basename=('WT_' if role=='wt_reference' else '')+f'input_1_{repeat-1}.pdb'
            require(path.as_posix()==(oldroot/'jobs'/oldjob/basename).as_posix(),'Source does not preserve corresponding original repeat WT reference')
            require(path.as_posix() in protocol['source_hashes'],'Source lacks frozen hash')
            validate_source((root/path).read_text(encoding='ascii'),sequence)
            source_uid=f'{uid}_{role}'
            row=dict(source_uid=source_uid,pair_uid=uid,old_job_id=oldjob,peptide_id=pid,parent_id='P53_WT',seed=seed,model_number=model,repeat=repeat,role=role,sequence=sequence,experimental_direction=direction,source_path=path.as_posix(),source_sha256=protocol['source_hashes'][path.as_posix()],input_path=f'sources/{source_uid}/input.pdb')
            sources.append(row)
            for condition in ('fixed','flexible'): jobs.append(row|dict(job_id=f'{source_uid}_{condition}',condition=condition))
    require(len(sources)==120 and len(jobs)==240,'Invalid source/job manifest count')
    return sources,jobs


def prepare(output=OUTPUT,root=ROOT) -> dict:
    output,root=Path(output).resolve(),Path(root).resolve()
    protocol=load_json(output/'protocol.json'); validate_protocol(protocol)
    check_hashes(protocol['source_hashes'],root); verify_environment(protocol)
    sources,jobs=build_manifests(protocol,root)
    prepared=output/'preparation.json'
    if prepared.exists():
        p=load_json(prepared)
        require(p['protocol_sha256']==sha(output/'protocol.json'),'Frozen protocol changed')
        require(p['implementation_hashes']==implementation_hashes(root),'Frozen implementation changed')
        check_hashes(p['prepared_hashes'],output)
        require((output/'source_manifest.csv').read_text()==csv_text(sources),'Source manifest identity changed')
        require((output/'job_manifest.csv').read_text()==csv_text(jobs),'Job manifest identity changed')
        return p
    require(not any((output/x).exists() for x in ('sources','jobs','source_manifest.csv','job_manifest.csv')),'Untracked preparation files; refusing overwrite')
    hashes={}
    for row in sources:
        destination=output/row['input_path']; destination.parent.mkdir(parents=True)
        shutil.copyfile(root/row['source_path'],destination)
        require(sha(destination)==row['source_sha256'],'Source copy changed')
        hashes[row['input_path']]=sha(destination)
    for filename,rows in [('source_manifest.csv',sources),('job_manifest.csv',jobs)]:
        (output/filename).write_text(csv_text(rows),encoding='utf-8',newline='\n'); hashes[filename]=sha(output/filename)
    p=dict(schema_version=1,protocol_sha256=sha(output/'protocol.json'),implementation_hashes=implementation_hashes(root),prepared_hashes=hashes,source_count=120,structure_job_count=240,pair_count=120,contrast_count=60)
    save_json(prepared,p)
    save_json(output/'execution_state.json',dict(schema_version=1,preparation_sha256=sha(prepared),jobs={},hydrations={},invalid_sources={}))
    return p


def artifacts(job: Path) -> dict:
    return {p.relative_to(job).as_posix():sha(p) for p in sorted(job.rglob('*')) if p.is_file() and p.name!='result.json' and 'attempts' not in p.relative_to(job).parts}


def archive_attempt(job: Path):
    paths=[x for x in job.iterdir() if x.name!='attempts'] if job.exists() else []
    if not paths: return None
    archive=job/'attempts'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    archive.mkdir(parents=True)
    for path in paths: shutil.move(str(path),str(archive/path.name))
    save_json(archive/'artifact_hashes.json',{p.relative_to(archive).as_posix():sha(p) for p in archive.rglob('*') if p.is_file()})
    return archive


def archive_incomplete_hydration(source: Path):
    common=source/'common'
    require(common.is_dir() and not (common/'hydration.json').exists(),'Completed/untracked hydration cannot be regenerated')
    archive=source/'hydration_attempts'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    archive.parent.mkdir(parents=True,exist_ok=True)
    shutil.move(str(common),str(archive))
    save_json(archive/'artifact_hashes.json',{p.relative_to(archive).as_posix():sha(p) for p in archive.rglob('*') if p.is_file()})
    return archive


def recover_untracked_hydration(source: Path,common: Path,row: dict) -> str:
    """Accept a completed worker cutpoint only after independent identity checks.

    The worker may finish its atomic hydration metadata write immediately before
    the parent is interrupted. Missing metadata is an incomplete attempt; changed
    completed metadata/artifacts must never trigger silent regeneration.
    """
    from pipeline.backbone_relax import strip_hydrogens, mapping_keys
    import numpy as np
    require(sha(source)==row['source_sha256'],'Hydration recovery copied source hash changed')
    metadata=load_json(common/'hydration.json')
    require(metadata['source_uid']==row['source_uid'] and metadata['source_sha256']==row['source_sha256'],'Completed hydration source identity/hash changed')
    expected_artifacts={'heavy_start.pdb','hydrated.pdb','positions.json'}
    require(set(metadata['artifact_hashes'])==expected_artifacts,'Completed hydration artifact manifest is incomplete/ambiguous')
    require({p.name for p in common.iterdir() if p.is_file()}==expected_artifacts|{'hydration.json'},'Untracked completed hydration artifacts')
    check_hashes(metadata['artifact_hashes'],common)
    require(metadata['pH']==7.4 and metadata['platform']=='CPU' and metadata['platform_properties']=={'Threads':'1','DeterministicForces':'true'},'Completed hydration settings changed')
    digest=hashlib.sha256(f"20261007:{row['source_uid']}:{row['source_sha256']}".encode()).digest()
    require(metadata['hydrogen_seed']==int.from_bytes(digest[:8],'big')%(2**32),'Completed hydration deterministic seed changed')
    require(math.isfinite(metadata['heavy_max_displacement_A']) and metadata['heavy_max_displacement_A']<=1e-9,'Completed hydration moved heavy atoms')
    source_text=source.read_text(encoding='ascii')
    heavy_text=(common/'heavy_start.pdb').read_text(encoding='ascii')
    require(strip_hydrogens(source_text)==heavy_text,'Completed hydration heavy starting template changed')
    original=validate_source(heavy_text,row['sequence'])
    hydrated_text=(common/'hydrated.pdb').read_text(encoding='ascii')
    validate_source(hydrated_text,row['sequence'])
    data=load_json(common/'positions.json'); atoms=data['atoms']; keys=mapping_keys(atoms)
    positions=np.asarray(data['positions'],dtype=float)
    require(data['units']=='angstrom' and positions.shape==(len(atoms),3) and np.isfinite(positions).all(),'Completed hydration full precision coordinates invalid')
    require(len(keys)==len(set(keys)) and [a['index'] for a in atoms]==list(range(len(atoms))),'Completed hydration atom mapping ambiguous')
    heavy={key:tuple(positions[i]) for i,key in enumerate(keys) if atoms[i]['element'] not in ('H','D')}
    require(geo.displacement(original,heavy)['max_A']<=1e-9,'Completed hydration full precision heavy coordinates changed')
    pdb_atoms=[line for line in hydrated_text.splitlines() if line.startswith('ATOM  ')]
    require(len(pdb_atoms)==len(atoms),'Completed hydration PDB/full precision atom count differs')
    for i,(line,atom) in enumerate(zip(pdb_atoms,atoms)):
        require((line[21],int(line[22:26]),line[12:16].strip())==keys[i] and line[17:20].strip()==atom['residue_name'] and line[76:78].strip().upper()==atom['element'].upper(),'Completed hydration PDB atom mapping differs')
        xyz=np.asarray([float(line[k:k+8]) for k in (30,38,46)])
        require(np.linalg.norm(xyz-positions[i])<=.000867,'Completed hydration PDB coordinates differ beyond rounding')
    return sha(common/'hydration.json')


def invalidate_job(job: Path,job_id: str,signature: str,reason: str) -> dict:
    job.mkdir(parents=True,exist_ok=True); archive_attempt(job)
    (job/'permanently_invalid.txt').write_text(reason+'\n',encoding='utf-8')
    result=dict(job_id=job_id,signature=signature,status='failed',provenance_valid=False,failure_reason=reason,artifact_hashes=artifacts(job))
    save_json(job/'result.json',result)
    return result


def require_promotable(job: Path):
    require(not (job/'permanently_invalid.txt').exists(),'Changed provenance is permanently invalid; cannot promote/retry')
    if (job/'result.json').exists(): require(load_json(job/'result.json').get('provenance_valid',False),'Previously invalid job cannot promote')


def signature(output: Path,row: dict):
    return hashlib.sha256((sha(output/'preparation.json')+'\n'+json.dumps(row,sort_keys=True)).encode()).hexdigest()


def validate_cached(job: Path,job_id: str,expected_signature: str,expected_hash: str):
    require(sha(job/'result.json')==expected_hash,'Tracked result hash changed')
    result=load_json(job/'result.json')
    require(result['job_id']==job_id and result['signature']==expected_signature,'Cached job identity/provenance changed')
    require(result['artifact_hashes']==artifacts(job),'Cached output artifact hashes changed')
    require(result['status'] in ('ok','failed') and type(result['provenance_valid']) is bool,'Invalid cached job state')
    if not result['provenance_valid']: require(result['status']=='failed' and 'raw' not in result,'Provenance failure cannot carry score')
    if result['status']=='failed': require(result['failure_reason'] and 'raw' not in result,'Failed job cannot carry score'); return result
    require(result['provenance_valid'],'Valid score with invalid provenance')
    require(load_json(job/'analyse.command.json')['returncode']==0,'Actual AnalyseComplex command failed')
    require(result['raw']==parse_score(job,'relaxed.pdb'),'Cached score differs from raw tables')
    require(result['relaxation']==load_json(job/'relaxation.json'),'Relaxation metadata differs from raw artifact')
    return result


def verify_coordinates(output: Path,row: dict,result: dict):
    from pipeline.backbone_relax import mapping_keys, frozen_indices
    common=output/'sources'/row['source_uid']/'common'; job=output/'jobs'/row['job_id']
    metadata=load_json(common/'hydration.json'); check_hashes(metadata['artifact_hashes'],common)
    require(metadata['source_uid']==row['source_uid'] and metadata['source_sha256']==row['source_sha256'],'Common hydration identity changed')
    relaxation=result['relaxation']
    require(relaxation['hydration_sha256']==sha(common/'hydration.json') and relaxation['start_positions_sha256']==sha(common/'positions.json'),'Condition does not share tracked common hydration')
    initial=load_json(common/'positions.json'); final=load_json(job/'positions.json')
    require(initial['atoms']==final['atoms'] and final['start_positions_sha256']==sha(common/'positions.json'),'Full precision atom mapping/start mismatch')
    import numpy as np
    start=np.asarray(initial['positions']); stop=np.asarray(final['positions'])
    require(np.isfinite(start).all() and np.isfinite(stop).all() and start.shape==stop.shape,'Invalid full precision positions')
    indices=frozen_indices(mapping_keys(initial['atoms']),row['condition'])
    maximum=float(np.linalg.norm(stop[indices]-start[indices],axis=1).max())
    require(maximum<=1e-6 and math.isclose(maximum,relaxation['frozen_max_A'],abs_tol=1e-12),'Actual full precision frozen atom movement')
    original=validate_source((output/row['input_path']).read_text(),row['sequence'])
    heavy=validate_source((job/'relaxed.pdb').read_text(),row['sequence'])
    require(original.keys()==heavy.keys(),'Heavy atom identities changed during relaxation/scoring')
    stop_heavy={k:tuple(stop[i]) for i,k in enumerate(mapping_keys(final['atoms'])) if final['atoms'][i]['element'] not in ('H','D')}
    require(geo.displacement(stop_heavy,heavy)['max_A']<=.000867,'Scored PDB differs from full precision result beyond rounding')


def worker_hydrate(output: Path,root: Path,source_uid: str):
    from pipeline.backbone_relax import hydrate
    rows=read_csv(output/'source_manifest.csv'); row=next(x for x in rows if x['source_uid']==source_uid)
    p=load_json(output/'preparation.json'); require(p['implementation_hashes']==implementation_hashes(root),'Worker implementation changed')
    require(sha(output/row['input_path'])==row['source_sha256'],'Hydration input changed')
    hydrate(output/row['input_path'],output/'sources'/source_uid/'common',source_uid,row['sequence'])


def worker_job(output: Path,root: Path,job_id: str):
    from pipeline.backbone_relax import relax
    row=next(x for x in read_csv(output/'job_manifest.csv') if x['job_id']==job_id)
    job=output/'jobs'/job_id; job.mkdir(parents=True,exist_ok=True); require_promotable(job)
    sig=signature(output,row)
    result=dict(job_id=job_id,source_uid=row['source_uid'],signature=sig,status='failed',provenance_valid=True,failure_reason='Worker did not complete')
    try:
        protocol=load_json(output/'protocol.json'); p=load_json(output/'preparation.json')
        require(p['implementation_hashes']==implementation_hashes(root),'Worker implementation changed')
        check_hashes(protocol['source_hashes'],root); verify_environment(protocol)
        require(sha(output/row['input_path'])==row['source_sha256'],'Copied source changed')
        result['relaxation']=relax(output/'sources'/row['source_uid']/'common',job,row['condition'],row['sequence'])
        foldx=root/'tools/foldx/foldx5_Linux/foldx_20261231'
        molecules=root/'tools/foldx/foldx5_Linux/molecules'
        targets=sorted(molecules.glob('*.json')); require(len(targets)==15,'Expected official 15 molecule JSONs')
        (job/'molecules').mkdir()
        for source in targets: shutil.copyfile(source,job/'molecules'/source.name)
        for source in targets:
            relative=source.relative_to(root).as_posix()
            require(sha(job/'molecules'/source.name)==protocol['source_hashes'][relative],'Frozen copied molecule resource changed')
        (job/'analyse.cfg').write_text(ANALYSE_CONFIG,encoding='ascii')
        scored_pdb_sha256=sha(job/'relaxed.pdb')
        command=analyse_command(foldx)
        log=command_log(command,job,'analyse',timeout=3600)
        require(log['returncode']==0,'AnalyseComplex subprocess failed: '+log.get('error',''))
        require(sha(job/'relaxed.pdb')==scored_pdb_sha256,'Scored PDB changed during AnalyseComplex')
        for source in targets:
            require(sha(job/'molecules'/source.name)==protocol['source_hashes'][source.relative_to(root).as_posix()],'Frozen molecule resource changed during AnalyseComplex')
        result['raw']=parse_score(job,'relaxed.pdb')
        result['scored_pdb_sha256']=scored_pdb_sha256
        result.update(status='ok',failure_reason='')
        verify_coordinates(output,row,result)
        check_hashes(protocol['source_hashes'],root); verify_environment(protocol)
    except Exception as error:
        result.update(status='failed',failure_reason=f'{type(error).__name__}: {error}')
        result.pop('raw',None)
        if any(word in str(error).lower() for word in ('changed','provenance','hash','identity','mapping','frozen atom')):
            result['provenance_valid']=False
            (job/'permanently_invalid.txt').write_text(result['failure_reason']+'\n')
    result['artifact_hashes']=artifacts(job); save_json(job/'result.json',result)
    return result


def worker_command(output: Path,root: Path,kind: str,uid: str):
    return [sys.executable,str(root/'pipeline/backbone_diagnostic.py'),kind,'--output',str(output),'--root',str(root),'--uid',uid]


def collect(output: Path,jobs: list[dict],mutate=False):
    state=load_json(output/'execution_state.json'); results={}
    require(state['preparation_sha256']==sha(output/'preparation.json'),'Execution preparation hash changed')
    for row in jobs:
        uid=row['job_id']; job=output/'jobs'/uid
        if uid not in state['jobs']: continue
        try:
            result=validate_cached(job,uid,signature(output,row),state['jobs'][uid]['result_sha256'])
            if result['status']=='ok': verify_coordinates(output,row,result)
        except Exception as error:
            if not mutate: raise
            result=invalidate_job(job,uid,signature(output,row),'Changed provenance: '+str(error))
            state['jobs'][uid]={'result_sha256':sha(job/'result.json')}
        results[uid]=result
    if mutate: save_json(output/'execution_state.json',state)
    return results


def run(output=OUTPUT,root=ROOT,phase='pilot',workers=2):
    output,root=Path(output).resolve(),Path(root).resolve(); require(workers in (1,2),'At most two CPU workers')
    prepare(output,root)
    alljobs=read_csv(output/'job_manifest.csv'); selected=[x for x in alljobs if phase=='full' or int(x['repeat'])==1]
    existing=collect(output,alljobs,mutate=True)
    if phase=='full':
        pilot_rows,_,_=analysis.assemble_pairs([x for x in alljobs if int(x['repeat'])==1],existing)
        require(analysis.technical_gate(pilot_rows,48)['technical_pass'],'Full run requires hash-verified technical pilot pass; biological score direction is not a gate')
    state=load_json(output/'execution_state.json')
    source_rows={x['source_uid']:x for x in selected}
    for source_uid in sorted({x['source_uid'] for x in selected}):
        source=output/'sources'/source_uid; common=source/'common'
        if source_uid in state['invalid_sources']: continue
        try:
            if source_uid in state['hydrations']:
                require(sha(common/'hydration.json')==state['hydrations'][source_uid],'Tracked common hydration changed')
                recover_untracked_hydration(source/'input.pdb',common,source_rows[source_uid])
            else:
                if (common/'hydration.json').exists():
                    state['hydrations'][source_uid]=recover_untracked_hydration(source/'input.pdb',common,source_rows[source_uid])
                else:
                    if common.exists(): archive_incomplete_hydration(source)
                    log=command_log(worker_command(output,root,'_hydrate',source_uid),source,'hydrate',timeout=3600)
                    require(log['returncode']==0 and (common/'hydration.json').exists(),'Common hydration process failed')
                    state['hydrations'][source_uid]=recover_untracked_hydration(source/'input.pdb',common,source_rows[source_uid])
        except Exception as error:
            state['invalid_sources'][source_uid]=str(error)
        save_json(output/'execution_state.json',state)
    pending=[]
    for row in selected:
        uid=row['job_id']; job=output/'jobs'/uid
        if uid in existing: continue
        if row['source_uid'] in state['invalid_sources']:
            invalidate_job(job,uid,signature(output,row),'Source hydration/provenance failure: '+state['invalid_sources'][row['source_uid']]); state['jobs'][uid]={'result_sha256':sha(job/'result.json')}; continue
        job.mkdir(parents=True,exist_ok=True); require_promotable(job); archive_attempt(job); pending.append(row)
    save_json(output/'execution_state.json',state)
    def process(row):
        job=output/'jobs'/row['job_id']; start=time.monotonic()
        # Launch logs are outside the worker's artifact directory, avoiding write/hash races.
        command=worker_command(output,root,'_worker',row['job_id'])
        proc=subprocess.run(command,capture_output=True,text=True,encoding='utf-8',errors='replace',shell=False)
        launches=output/'launches'; launches.mkdir(exist_ok=True)
        save_json(launches/(row['job_id']+'.json'),dict(command=command,returncode=proc.returncode,stdout=proc.stdout,stderr=proc.stderr,seconds=time.monotonic()-start))
        if proc.returncode!=0 or not (job/'result.json').exists():
            result=dict(job_id=row['job_id'],source_uid=row['source_uid'],signature=signature(output,row),status='failed',provenance_valid=True,failure_reason=f'Worker process failed ({proc.returncode}): {proc.stderr[-2000:]}',artifact_hashes=artifacts(job)); save_json(job/'result.json',result)
        return row['job_id']
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for uid in pool.map(process,pending):
            state['jobs'][uid]={'result_sha256':sha(output/'jobs'/uid/'result.json')}; save_json(output/'execution_state.json',state)
            print(f'{phase}: {len(state["jobs"])} tracked structure jobs',flush=True)
    result=analysis.analyze(output,alljobs,collect(output,alljobs,mutate=True),phase)
    return result


def verify(output=OUTPUT,root=ROOT):
    output,root=Path(output).resolve(),Path(root).resolve(); prepare(output,root)
    jobs=read_csv(output/'job_manifest.csv'); results=collect(output,jobs)
    state=load_json(output/'execution_state.json')
    for uid,expected in state['hydrations'].items():
        common=output/'sources'/uid/'common'; require(sha(common/'hydration.json')==expected,'Common hydration hash changed'); check_hashes(load_json(common/'hydration.json')['artifact_hashes'],common)
    published=output/'analysis/published_artifact_hashes.json'
    if published.exists(): check_hashes(load_json(published),output)
    rows,pairs,_=analysis.assemble_pairs(jobs,results)
    gate=analysis.technical_gate(rows,240)
    gate.update(verification='all tracked actual outputs/hash/coordinate checks passed',complete=gate['technical_pass'] and len(pairs)==120 and all(x['status']=='ok' for x in pairs),tracked_jobs=len(results))
    return gate


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('prepare','baseline','run','analyze','verify','_worker','_hydrate'))
    parser.add_argument('--output',type=Path,default=OUTPUT); parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--phase',choices=('pilot','full'),default='full'); parser.add_argument('--workers',type=int,choices=(1,2),default=2); parser.add_argument('--uid')
    args=parser.parse_args(); output=args.output.resolve(); root=args.root.resolve()
    if args.command=='prepare': result=prepare(output,root)
    elif args.command=='run': result=run(output,root,args.phase,args.workers)
    elif args.command=='_worker': result=worker_job(output,root,args.uid)
    elif args.command=='_hydrate': worker_hydrate(output,root,args.uid); result={'hydrated':args.uid}
    elif args.command=='baseline': prepare(output,root); result=analysis.baseline(output,root,load_json(output/'protocol.json'))
    elif args.command=='analyze': prepare(output,root); jobs=read_csv(output/'job_manifest.csv'); result=analysis.analyze(output,jobs,collect(output,jobs,mutate=True),args.phase)
    else: result=verify(output,root)
    print(json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False))


if __name__=='__main__': main()
