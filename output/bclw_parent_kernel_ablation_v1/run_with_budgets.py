from __future__ import annotations
import hashlib, importlib.metadata as md, json, os, pathlib, shutil, subprocess, sys, time
from datetime import datetime, timezone

WORK=pathlib.Path('/mnt/e/Users/graphpep与pepmlm')
OUT=WORK/'output/bclw_parent_kernel_ablation_v1'
AUDIT=WORK/'output/structure_input_audit_bclw_20261020'
PILOT=WORK/'output/bh3_bclw_six_complex_pilot_20261020'
PLAN_PATH=AUDIT/'next_compute_protocol_prepared_v1.json'
PLAN=json.loads(PLAN_PATH.read_text())
LIMIT_GPU=PLAN['budget']['gpu_tagged_process_wall_seconds_max']
LIMIT_CPU=PLAN['budget']['cpu_relax_wall_seconds_max']
LIMIT_DISK=PLAN['budget']['disk_bytes_max']
EVENTS=[]
GPU_SECONDS=0.0
CPU_SECONDS=0.0


def now(): return datetime.now(timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def tree_bytes():
    total=0
    for p in OUT.rglob('*'):
        if p.is_file():
            try: total+=p.stat().st_size
            except FileNotFoundError: pass
    return total
def write_json(name,obj):
    p=OUT/name; t=p.with_suffix(p.suffix+'.tmp')
    t.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    t.replace(p)
def event(name,status,**extra):
    EVENTS.append({'time_utc':now(),'name':name,'status':status,**extra})
    write_json('execution_record.json',{'protocol_id':PLAN['protocol_id'],'events':EVENTS,
                'gpu_prediction_wall_seconds':GPU_SECONDS,'cpu_relax_wall_seconds':CPU_SECONDS,
                'disk_bytes_current':tree_bytes(),'disk_bytes_max':LIMIT_DISK})

def finite_tree(x):
    if isinstance(x,float): return x==x and abs(x)!=float('inf')
    if isinstance(x,dict): return all(finite_tree(v) for v in x.values())
    if isinstance(x,list): return all(finite_tree(v) for v in x)
    return True

def check_originals():
    manifest=json.loads((PILOT/'run_result_manifest_v1.json').read_text())
    mismatches=[]
    for row in manifest['files']:
        p=PILOT/row['path']
        if not p.is_file() or p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:
            mismatches.append(row['path'])
    return {'count':len(manifest['files']),'all_match':len(manifest['files'])==265 and not mismatches,'mismatches':mismatches}

def preflight():
    s=PLAN['constant_prediction_settings']
    checks={}
    checks['output_root_clean_of_results']=OUT.exists() and set(x.name for x in OUT.iterdir())=={'run_with_budgets.py','logs','prepared_protocol.json'} and not any((OUT/'logs').iterdir())
    checks['input_sha256']=sha(PLAN['input_fasta'])
    checks['input_matches_frozen']=checks['input_sha256']==PLAN['input_fasta_sha256']
    checks['batch_source_sha256']=sha('/home/g/miniconda3/envs/colabfold_exp0/lib/python3.11/site-packages/colabfold/batch.py')
    checks['weights_sha256']=sha(s['model_parameter_file'])
    checks['relax_script_sha256']=sha(WORK/'pipeline/relax_feedback_panel.py')
    checks['pilot_manifest']=check_originals()
    checks['versions']={n:md.version(n) for n in ('colabfold','jax','jaxlib','openmm')}
    checks['gpu']={}
    for key,cmd in {
      'gpu_snapshot':['nvidia-smi'],
      'gpu_memory':['nvidia-smi','--query-gpu=name,memory.total,memory.free,memory.used,utilization.gpu','--format=csv,noheader'],
      'gpu_compute_processes':['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],
      'existing_target_processes':['bash','-lc','pgrep -af "colabfold_batch|relax_feedback_panel.py" || true']}.items():
        r=subprocess.run(cmd,capture_output=True,text=True,timeout=30)
        checks['gpu'][key]={'returncode':r.returncode,'stdout':r.stdout.strip(),'stderr':r.stderr.strip()}
    checks['existing_target_processes']=[]
    for procdir in pathlib.Path('/proc').iterdir():
        if not procdir.name.isdigit() or int(procdir.name) in (os.getpid(),os.getppid()): continue
        try: cmd=(procdir/'cmdline').read_bytes().replace(b'\\x00',b' ').decode(errors='replace')
        except (FileNotFoundError,PermissionError,ProcessLookupError): continue
        if 'colabfold_batch' in cmd or 'relax_feedback_panel.py' in cmd: checks['existing_target_processes'].append({'pid':int(procdir.name),'cmdline':cmd[:500]})
    st=os.statvfs('/mnt/e')
    checks['disk_free_bytes']=st.f_bavail*st.f_frsize
    checks['prepared_protocol_sha256']=sha(PLAN_PATH)
    if (checks['batch_source_sha256']!='a65a1fb4301ae82de8b307a8d698f327a79b8fa66c8f053f4ddeff52c87020c4' or
        checks['weights_sha256']!=s['model_parameter_sha256'] or
        checks['relax_script_sha256']!='f88c3489bf89eb279a548f500eaa4102d04f0a0113bc73271f7e74c8d07703d8' or
        not checks['input_matches_frozen'] or not checks['pilot_manifest']['all_match'] or
        checks['versions']!=s['versions'] or not checks['output_root_clean_of_results'] or
        checks['disk_free_bytes']<LIMIT_DISK or
        checks['gpu']['gpu_compute_processes']['stdout'] or checks['existing_target_processes']):
        raise RuntimeError('Preflight mismatch; no prediction launched: '+json.dumps(checks,ensure_ascii=False))
    write_json('preflight.json',{'time_utc':now(),'protocol_id':PLAN['protocol_id'],'checks':checks,
                                'budgets':PLAN['budget'],'frozen_quality_gate':PLAN['quality_gate_unchanged'],
                                'additional_diagnostic_gate':PLAN['additional_diagnostic_gate']})
    event('preflight','passed',input_sha256=checks['input_sha256'],pilot_files=checks['pilot_manifest']['count'],gpu_memory=checks['gpu']['gpu_memory']['stdout'],disk_free_bytes=checks['disk_free_bytes'])

def run_limited(label,argv,env_extra,category,limit):
    global GPU_SECONDS, CPU_SECONDS
    if tree_bytes()>LIMIT_DISK: raise RuntimeError('Disk cap already exceeded before '+label)
    if category=='gpu' and GPU_SECONDS>=LIMIT_GPU: raise RuntimeError('GPU wall budget exhausted before '+label)
    if category=='cpu' and CPU_SECONDS>=LIMIT_CPU: raise RuntimeError('CPU wall budget exhausted before '+label)
    env=os.environ.copy(); env.update(env_extra)
    log=OUT/'logs'/f'{label}.log'; log.parent.mkdir(parents=True,exist_ok=True)
    start= time.monotonic(); pid=None
    event(label,'started',argv=argv,environment=env_extra,wall_limit_remaining_seconds=limit-(GPU_SECONDS if category=='gpu' else CPU_SECONDS))
    with log.open('wb') as stream:
        proc=subprocess.Popen(argv,cwd=str(WORK),env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        pid=proc.pid
        last=0
        while proc.poll() is None:
            elapsed=time.monotonic()-start
            if category=='gpu': GPU_SECONDS=elapsed+sum(e.get('wall_seconds',0) for e in EVENTS if e.get('category')=='gpu' and e.get('name')!=label)
            else: CPU_SECONDS=elapsed+sum(e.get('wall_seconds',0) for e in EVENTS if e.get('category')=='cpu' and e.get('name')!=label)
            size=tree_bytes()
            if elapsed>limit or size>LIMIT_DISK:
                proc.terminate()
                try: proc.wait(timeout=10)
                except subprocess.TimeoutExpired: proc.kill(); proc.wait()
                spent=time.monotonic()-start
                if category=='gpu': GPU_SECONDS+=spent
                else: CPU_SECONDS+=spent
                event(label,'stopped_budget',category=category,wall_seconds=spent,disk_bytes=size,pid=pid)
                raise RuntimeError(f'{label} stopped at a frozen budget cap')
            if int(elapsed//30)>last:
                last=int(elapsed//30); print(f'[{label}] elapsed={int(elapsed)}s disk={size}B',flush=True)
            time.sleep(1)
        spent=time.monotonic()-start
    if category=='gpu': GPU_SECONDS+=spent
    else: CPU_SECONDS+=spent
    if tree_bytes()>LIMIT_DISK: event(label,'stopped_budget',category=category,wall_seconds=spent,pid=pid); raise RuntimeError('disk cap exceeded')
    if proc.returncode!=0:
        event(label,'failed',category=category,returncode=proc.returncode,wall_seconds=spent,pid=pid,log=str(log))
        raise RuntimeError(f'{label} failed exit={proc.returncode}; see {log}')
    event(label,'completed',category=category,returncode=0,wall_seconds=spent,pid=pid,log=str(log))
    return spent

def validate_prediction(seed,command):
    out=OUT/f'seed_{seed}'/'structures'
    pfiles=sorted(out.glob('*_unrelaxed_rank_*.pdb'))
    scores=sorted(out.glob('*_scores_rank_*.json'))
    configs=out/'config.json'
    if len(pfiles)!=1 or len(scores)!=1 or not configs.is_file(): raise RuntimeError(f'seed {seed}: expected one PDB, one scores JSON, and config; got {len(pfiles)}, {len(scores)}, {configs.exists()}')
    def reject_constant(v): raise ValueError('non-finite JSON constant '+v)
    score=json.loads(scores[0].read_text(),parse_constant=reject_constant)
    if not finite_tree(score): raise RuntimeError(f'seed {seed}: non-finite score JSON')
    cfg=json.loads(configs.read_text(),parse_constant=reject_constant)
    baseline=json.loads((PILOT/f'seed_{seed}'/'structures'/'config.json').read_text())
    changed={k:{'baseline':baseline.get(k,'<absent>'),'new':cfg.get(k,'<absent>')} for k in sorted(set(baseline)|set(cfg)) if baseline.get(k,'<absent>')!=cfg.get(k,'<absent>')}
    expected={'msa_mode':'single_sequence','use_templates':False,'num_relax':0,'model_type':'alphafold2_multimer_v3','num_models':1,'num_recycles':3,'recycle_early_stop_tolerance':0.0,'num_ensemble':1,'model_order':[1],'pair_mode':'unpaired','random_seed':seed,'num_seeds':1,'use_fast_kernels':False,'max_seq':1}
    failed={k:{'expected':v,'actual':cfg.get(k,'<absent>')} for k,v in expected.items() if cfg.get(k)!=v}
    if failed: raise RuntimeError(f'seed {seed}: effective config mismatch: '+json.dumps(failed))
    # All changes to same-seed effective config are recorded; only query count (panel to parent-only) and kernel flag are expected.
    unexpected=sorted(set(changed)-{'num_queries','use_fast_kernels'})
    event(f'config_check_seed_{seed}','passed' if not failed else 'failed',changed_fields=changed,unexpected_changed_fields=unexpected)
    return {'pdb':str(pfiles[0]),'pdb_sha256':sha(pfiles[0]),'scores':str(scores[0]),'scores_sha256':sha(scores[0]),'config':str(configs),'config_sha256':sha(configs),'changed_fields':changed,'unexpected_changed_fields':unexpected}

def main():
    if OUT.exists():
        if set(p.name for p in OUT.iterdir())!={'run_with_budgets.py','logs','prepared_protocol.json'} or any((OUT/'logs').iterdir()): raise RuntimeError('Output root contains unexpected files at launch; refusing overwrite')
    else: OUT.mkdir(parents=True)
    (OUT/'logs').mkdir(exist_ok=True)
    shutil.copy2(PLAN_PATH,OUT/'prepared_protocol.json')
    preflight()
    results=[]
    try:
      for cmd in PLAN['commands']:
        seed=cmd['seed']
        if sha(PLAN['input_fasta'])!=PLAN['input_fasta_sha256']: raise RuntimeError('Input changed between seeds')
        run_limited(f'predict_seed_{seed}',cmd['argv'],cmd['environment'],'gpu',LIMIT_GPU-GPU_SECONDS)
        pred=validate_prediction(seed,cmd)
        # Catch effective derived settings, but do not alter or rerun the experimental branch.
        results.append({'seed':seed,'prediction':pred})
        run_limited(f'relax_seed_{seed}',cmd['relax_argv'],{'JAX_PLATFORMS':'cpu','XLA_PYTHON_CLIENT_PREALLOCATE':'false','OPENMM_CPU_THREADS':'2'},'cpu',LIMIT_CPU-CPU_SECONDS)
        relax_manifest=OUT/f'seed_{seed}'/'relaxed'/'relaxation_manifest.json'
        if not relax_manifest.is_file(): raise RuntimeError(f'seed {seed}: relaxation manifest missing')
        rm=json.loads(relax_manifest.read_text())
        if len(rm.get('models',[]))!=1 or rm.get('models',[{}])[0].get('settings')!=PLAN['constant_relax_settings']:
            raise RuntimeError(f'seed {seed}: relaxation settings/row count differ from frozen protocol')
        results[-1]['relaxation_manifest']=str(relax_manifest)
        results[-1]['relaxation_manifest_sha256']=sha(relax_manifest)
        write_json('completed_jobs.json',{'protocol_id':PLAN['protocol_id'],'jobs':results,'gpu_prediction_wall_seconds':GPU_SECONDS,'cpu_relax_wall_seconds':CPU_SECONDS})
      event('all_jobs','completed',job_count=len(results),gpu_prediction_wall_seconds=GPU_SECONDS,cpu_relax_wall_seconds=CPU_SECONDS)
    except Exception as exc:
      event('run','failed',error_type=type(exc).__name__,error=str(exc),gpu_prediction_wall_seconds=GPU_SECONDS,cpu_relax_wall_seconds=CPU_SECONDS)
      write_json('failure.json',{'time_utc':now(),'error_type':type(exc).__name__,'error':str(exc),'gpu_prediction_wall_seconds':GPU_SECONDS,'cpu_relax_wall_seconds':CPU_SECONDS,'disk_bytes':tree_bytes()})
      raise
    finally:
      try:
        after=check_originals()
        write_json('postflight_original_integrity.json',after)
      except Exception as e: write_json('postflight_original_integrity.json',{'check_error':repr(e)})

if __name__=='__main__': main()