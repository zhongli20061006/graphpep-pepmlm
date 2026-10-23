from pathlib import Path
import csv, hashlib, json, random
from collections import Counter
import numpy as np
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
R=Path('/mnt/e/Users/graphpep与pepmlm'); O=R/'output/algorithm_hypothesis_plan_20261022'
meta_path=R/'output/feedback_panel_mdm2_20261005/parent_metadata.json'
meta=json.loads(meta_path.read_text()); target=meta['receptor_sequence']
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as stream:
  for b in iter(lambda:stream.read(1024*1024),b''): h.update(b)
 return h.hexdigest()
def wj(name,obj): (O/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
used=[meta_path,R/'pipeline/feedback_control_1ycr.py',R/'pipeline/step1_generate.py',R/'output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json',R/'output/reward_direction_source_screen_20261020/report_zh.md',R/'output/reward_direction_source_screen_20261020/raw/skempi_v2.csv',R/'output/bclw_parent_kernel_ablation_v1/decision.json']
rows=[]; parser=PDBParser(QUIET=True)
for p in meta['parents'][:2]:
 path=Path(p['native_path']); used.append(path); model=parser.get_structure(p['peptide_id'],str(path))[0]
 residues=[r for r in model['B'] if r.id[0]==' ']; observed=''.join(seq1(r.resname) for r in residues)
 full=p['sequence']; assert full.startswith(observed)
 a=[a.coord for r in model['A'] if r.id[0]==' ' for a in r if a.element.upper() not in ('H','D')]
 b=[a.coord for r in residues for a in r if a.element.upper() not in ('H','D')]
 bi=[i+1 for i,r in enumerate(residues) for a in r if a.element.upper() not in ('H','D')]
 dist2=np.sum((np.asarray(a)[:,None,:]-np.asarray(b)[None,:,:])**2,axis=-1)
 hits=np.argwhere(dist2<=36.0)
 contacting=sorted({int(bi[j]) for i,j in hits}); protected=[3,7,10]
 protocol_exclusions=[11] if p['peptide_id']=='P53WT_1YCR' else []
 eligible=sorted(set(contacting)-set(protected)-set(protocol_exclusions))
 rows.append({'parent_id':p['peptide_id'],'full_sequence':full,'full_sequence_length':len(full),'observed_coordinate_sequence':observed,'observed_coordinate_length':len(observed),'unresolved_positions':list(range(len(observed)+1,len(full)+1)), 'native_pdb':str(path),'native_pdb_sha256':sha(path),'geometry_contact_positions_le6A':contacting,'protected_positions':protected,'protocol_excluded_positions':protocol_exclusions,'eligible_positions':eligible,'eligible_count':len(eligible),'receptor_sequence_for_PepMLM':target,'receptor_length':len(target),'role':'TECHNICAL_DEVELOPMENT_ONLY_NO_AFFINITY_LABELS','note':'Full parent sequence is not automatically identical to the historical 12mer assay. PMI terminal P is sequence-supported but unresolved and cannot be masked in this proposal.'})
controls=[{'id':'MDM2_target','sequence':target,'role':'target_context'}]
for seed in (2026102201,2026102202):
 values=list(target); random.Random(seed).shuffle(values); seq=''.join(values)
 assert seq!=target and Counter(seq)==Counter(target)
 controls.append({'id':f'shuffled_context_{seed}','sequence':seq,'role':'ARTIFICIAL_ORDER_CONTROL_NOT_EXPERIMENTAL_NONBINDER','seed':seed})
assert len({x['sequence'] for x in controls})==3
wj('development_input_inventory.json',{'scope':'Existing native PDB geometry plus source-supported sequences; no affinity labels are attached to these engineering inputs.','inputs':rows})
(O/'contexts_prepared.fasta').write_text(''.join(f">{x['id']}\n{x['sequence']}\n" for x in controls))
wj('context_controls.json',{'controls':controls,'claim_limit':'Composition- and length-matched shuffled sequences test context order sensitivity; they do not establish biological specificity and are not experimental nonbinding examples.'})
with (O/'development_input_inventory.csv').open('w',encoding='utf-8-sig',newline='') as f:
 fields=['parent_id','full_sequence','full_sequence_length','observed_coordinate_sequence','observed_coordinate_length','unresolved_positions','geometry_contact_positions_le6A','protected_positions','protocol_excluded_positions','eligible_positions','eligible_count','native_pdb','native_pdb_sha256','role']; writer=csv.DictWriter(f,fieldnames=fields); writer.writeheader()
 for row in rows: writer.writerow({k:(';'.join(map(str,row[k])) if isinstance(row[k],list) else row[k]) for k in fields})
model_dir=Path('/home/g/pepmlm_pipeline/models/PepMLM-650M'); assets=[]
for p in sorted(model_dir.iterdir()):
 if p.is_file() and p.suffix in ('.json','.txt','.bin','.safetensors'):
  assets.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}); print('hashed_asset',p.name,p.stat().st_size,flush=True)
wj('existing_model_asset_manifest.json',{'model_dir':str(model_dir),'runtime':'/home/g/miniconda3/envs/graphpep_env/bin/python','runtime_exists':Path('/home/g/miniconda3/envs/graphpep_env/bin/python').is_file(),'assets':assets,'weights_loaded':False})
wj('source_reuse_manifest.json',{'files':[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in used],'all_reads_only':True})
sites=sum(x['eligible_count'] for x in rows)
proposal={
 'proposal_id':'RECEPTOR_CONTEXT_MUTATION_PROPOSAL_V0', 'status':'DESIGN_PROPOSAL_INPUTS_PREPARED_NOT_IMPLEMENTED_NOT_RUN',
 'current_turn':{'new_model_inference':0,'new_structure_predictions':0,'new_reward_calls':0,'training':0,'new_knowledge_claim':'none'},
 'hypothesis':'For single peptide substitutions, a receptor-context contrast of PepMLM masked marginal log-odds may separate target-dependent preferences from context-insensitive sequence preferences. It is a candidate proposal signal, not calibrated affinity.',
 'why_distinct_from_current_pipeline':'Current guided arm selects positions by contact-count-weighted static GraphPep residue feedback, then samples PepMLM; the new proposal evaluates site-and-amino-acid jointly and includes a receptor-context ablation.',
 'formula':{'d_R_i_a':'log p_PepMLM(a | R, P masked only at i) - log p_PepMLM(P_i | R, P masked only at i)','h_i_a':'d_target_i_a - median(d_shuffled_1_i_a, d_shuffled_2_i_a)','null_controls':'Two fixed composition-preserving receptor permutations; artificial technical controls, not affinity negatives.','units':'natural-log units, not kcal/mol, Kd, or binding deltaG','no_label_fitting':True},
 'stage1_technical_contract':{
   'input_inventory':'development_input_inventory.json','contexts':'contexts_prepared.fasta','parents':len(rows),'eligible_mask_sites':sites,'context_count':3,
   'base_logical_masked_inputs':3*sites,'repeat_target_inputs_for_batch_order_check':sites,'total_logical_masked_inputs_max':4*sites,
   'per_mask_outputs':'One model distribution gives all 20 standard amino-acid logits. It is not 19 separate structure predictions.',
   'method_arms':['new context-contrast score H','unadjusted target-conditioned log-odds d','existing static GraphPep position ranking (reuse only)','uniform random ranking over the exact same eligible sites and amino-acid pool'],
   'historical_label_exposure':'All MDM2/PMI data are DEVELOPMENT. No fitting or threshold selection using known affinity directions; no rerun of P27A. Existing counterexamples are not a fresh test set.',
   'forbid_named_historical_mutant_promotions':['P53_P27A','P53_L22A','PMI_N8A','PMI_W7A'],
   'scope_of_comparison':'Correct indexing, numerical consistency, dependence on receptor order, control disagreement and candidate coverage; no affinity direction accuracy is calculated in stage 1.',
   'precision':'Proposed FP32 inference, eval mode, no gradients, TF32 disabled, batch size <=2. No precision fallback on OOM.',
   'numeric_comparison_tolerance_nats':0.0001,
   'suggestion_filter':'d_target, H_against_control_1 and H_against_control_2 must each exceed the numeric tolerance. Otherwise abstain; never resample controls to force recommendations.',
   'ranking':'Among eligible substitutions passing the filter, rank by H descending; one top suggestion per site, at most six distinct single-substitution suggestions per parent. No physical candidates or structures are generated in stage 1.',
   'cache_policy':'Do not use the old peptide-sequence-only PPL cache across contexts. Use fresh outputs in a new directory or a key containing checkpoint, receptor, peptide, masked position, precision and scoring definition.',
   'budget_proposed_not_authorized':{'gpu_process_wall_seconds_max':600,'disk_bytes_max':67108864,'new_structure_jobs_max':0,'GraphPep_calls_max':0,'FoldX_calls_max':0,'ESM2_calls_max':0,'training_steps':0},
   'stops':['input/checkpoint hash mismatch','nonfinite logits or bad mask indexing','WT-to-WT difference not zero within tolerance','batch/order differences exceed frozen tolerance','OOM or runtime error','wall/disk cap reached'],
   'interpretation':['If all receptor-context effects are at numerical tolerance, the proposed context signal is not demonstrated for these inputs; HOLD expansion.','If controls disagree, report it and abstain on the affected substitution.','If no suggestions remain, report zero coverage; do not change controls or thresholds to rescue.','Passing technical checks only means the prototype runs and has interpretable context sensitivity. It does not validate affinity, specificity, novelty or feedback superiority.']
 },
 'stage2_affinity_validation':{
   'status':'HOLD_LABEL_AND_PROTOCOL_GATES_NOT_MET',
   'baseline_policy':'Keep GraphPep, FoldX, ESM2 and contact features as peer candidate methods; compare H and unadjusted d on exactly the same eligible experimental labels. Preserve their frozen 0.05 tie bands.',
   'new_method_ties':'Numeric zero/tie rule must be frozen in the new protocol; it is not a biologically calibrated effect-size threshold.',
   'metrics':'Use frozen within-system/background-balanced direction accuracy, both class recalls, coverage and all failures/ties. Paired uncertainty for the method difference must be reported; no superiority claim if inconclusive.',
   'data_floor':'At least 10 eligible pairs, at least four in each direction and at least two parent backgrounds for a full single-system signal; engineering pilots do not replace this.',
   'falsifier':'A sufficiently qualified panel that shows H failing the frozen signal gate or failing to add useful evidence beyond unadjusted d does not support the context-correction hypothesis. Do not tune on its errors and call that same panel independent.',
   'structure_branch':'Any given-experimental-pose scoring has separate source/mapping/validation gates. Experimental PDB inputs do not require AFM ipTM/PAE; old AFM failures remain unchanged.',
   'data_work_budget':'At most one hour of targeted review of already identified missing inputs/labels. Reuse local SKEMPI and existing source reports. No duplicate SI downloads, broad database expansion, contact, payment or wet lab.'
 },
 'stage3_closed_loop':{'status':'HOLD_PENDING_SIGNAL_VALIDATION','future_design':'Same eligible site pool, one mutation per child, identical proposal/sequence budgets, matched random plus unadjusted-PepMLM ablation, fixed structural/reward evaluation. Re-evaluation by the same optimizing score alone is not proof of improved binding.'},
 'alternatives_considered':[{'approach':'Local physicochemical paired scorer','tradeoff':'Uses experimental poses directly, but current model failures cannot be repaired just by another matched-difference formula; additional modeling/labels needed.'},{'approach':'Context-contrast proposal using installed PepMLM','tradeoff':'Recommended first technical falsification: no training, no new weights or folding; artificial null contexts may be out of distribution and must be tested as an ablation.'},{'approach':'Train a new GNN or end-to-end generator','tradeoff':'Deferred: current independently qualified labels do not justify model capacity or evaluation claims.'}],
 'novelty':'Masked log-odds and contrastive ideas have prior literature. This is a project hypothesis/prototype proposal, not a verified novel algorithm.'}
wj('proposal.json',proposal)
print(json.dumps({'prepared_parents':[(x['parent_id'],x['eligible_positions']) for x in rows],'eligible_sites':sites,'future_logical_inputs_cap':4*sites,'inference_this_turn':0,'model_weights_loaded':False},ensure_ascii=False,indent=2))