"""Prepare tables and an unexecuted, bounded diagnostic protocol in this audit directory."""
from __future__ import annotations

import json
from pathlib import Path

from audit_inputs import OUT, PILOT, ROOT, SOURCE, sha, writecsv, writej


def main() -> None:
    rows = []
    def add(component: str, assay: str, model: str, status: str, evidence: str, limit: str) -> None:
        rows.append(dict(component=component, assay=assay, current_model=model, correspondence=status,
                         evidence=evidence, remaining_limit=limit))
    add('受体核心', '人Bcl-w 1-164；方法未给完整克隆序列', 'Q92843 canonical 1-164，164 aa',
        '名义区间一致；不是经克隆序列核验的完全相同分子',
        'Foight2014 main Methods; Dutta2010 PMC2896288 Methods; actual FASTA; UniProt Q92843',
        'Q133R为rs910332天然变体；旧结构多为R133，assay身份未知；不能擅自改成R133')
    add('受体标签', 'MBP经TEV切割后N端c-myc；完整接头/切割残基未恢复', '省略标签、接头和切割残基',
        '明确省略', 'Dutta2010 protein-expression methods; frozen input mapping', '不凭常见c-myc序列猜接头；没有标签导致低置信度的直接证据')
    add('亲本', 'Bim_A2eT荧光23-mer', 'RPEIWITQELRRIGDEFNAYYAR', '氨基酸核心一致',
        'Foight2014 Table1; frozen FASTA; 12 raw/relaxed PDB chain checks', '工程亲本不是天然WT；化学身份见端基行')
    add('增强候选', 'I3dL，工程亲本背景', 'I13L，23-mer，其余同亲本', '单点核心一致',
        'Foight2014 Tables1/2; FASTA/A3M/PDB', '源方向候选不等于本项目正式v2标签或reward验证')
    add('减弱候选', 'I2dM，工程亲本背景', 'I6M，23-mer，其余同亲本', '单点核心一致',
        'Foight2014 Tables1/2; FASTA/A3M/PDB', '同上；不跨亲本或assay合并')
    add('肽N端', 'N-5/6-fluorescein amidite', '标准氨基酸链；无染料原子', '明确省略',
        'Foight2014 Methods; input_representation_frozen_pilot_v1.json; PDB',
        '染料可能改变测量/结合；没有本组带/不带染料配对结构，不能定为预测失败根因')
    add('肽C端', 'C-terminal amide', '标准残基链；未显式建酰胺', '明确省略',
        'Foight2014 Methods; PDB; FoldX noCterm=false', '不是精确电荷/化学复现；AFM置信度不能从FoldX端基设置推断')
    add('实验条件', '25 C，pH7.8，50mM NaCl，EDTA，Triton，5%DMSO',
        'AFM未编码溶液条件；FoldX 298K/pH7.4/ionStrength0.15M', '代理设置差异',
        'Foight2014 Methods; frozen FoldX protocol', '不作完整ΔG(binding)，不把条件差异直接当AFM失败解释')
    add('聚集状态', 'Foight样品活性/寡聚比例未提供；Dutta2010制备曾报道混合状态', '1受体核心+1肽',
        '单体代理', 'Dutta2010 methods; actual A/B chains', 'Dutta2010比例不能移植成2014该批样品实测值')
    add('直接Kd vs竞争Ki', 'Table2为带荧光肽直接结合Kd', '作为外部来源方向候选', '只用直接Kd',
        'Foight2014 Table2 and Supplement Methods', '未标记N-acetyl/C-amide竞争Ki另存；不得混成同一标签')
    writecsv('assay_to_input_correspondence.csv', rows)

    structures = [
        ('1ZY3','NMR+分子对接；10模型','A:178 aa；缺当前Met1，P117V/Q133R；含尾/His标签','B:BID 20-mer；源注R203K','相关BH3结合参考','不是Bim或1-164精确assay；NMR/对接模型不能当精确晶体真值','https://www.rcsb.org/structure/1ZY3'),
        ('4CIM','X-ray 1.50 A','A/B:163 aa突变Bcl-w；C29S/H43G/M46L/A49I/A128E/Q133R；缺164','P/Q:21 aa突变Bcl-w自身BH3片段','排除为直接Bim pose来源','原文为蛋白水解短链及异常二聚构象；GraphPep训练结构暴露','https://pmc.ncbi.nlm.nih.gov/articles/PMC4276867/'),
        ('3FDL','X-ray 1.78 A','A:Bcl-xL工程/缺失构建，158 aa','B:天然Bim 26-mer；核心相对A2eT为T7A','跨受体同家族姿态参考','受体、长度、亲本和端基均不等于assay；不能转成Bcl-w实验证据','https://www.rcsb.org/structure/3FDL'),
        ('1MK3','solution NMR','A:178 aa；缺Met1，P117V/Q133R；标签/尾部','无独立肽链','受体构象参考','游离态/尾部占据与截尾assay不同；不能评价Bim复合物方向','https://www.rcsb.org/structure/1MK3'),
        ('1O0L','solution NMR','A:188 aa；含N端残基及尾；A128E/Q133R','无独立肽链','受体构象参考','尾部及突变不等于assay；不是复合物','https://www.rcsb.org/structure/1O0L'),
        ('2Y6W','X-ray','A/B:177 aa；含标签；C29S/A128E/Q133R','Bcl-w二聚体，无独立BH3肽','异常/交换构象参考','不能把另一个受体亚基当作Bim链','https://www.rcsb.org/structure/2Y6W'),
        ('4K5A','X-ray 1.50 A','entity1 A:184 aa；核心缺Met1，P117V/Q133R；SIFTS Q1RMX3/Bos taurus','entity2 B:DARPin D12','受体可结合构象参考','原论文称human但当前数据库映射bovine；标签/变体/伴侣不同；不是Bim复合物','https://www.rcsb.org/structure/4K5A'),
        ('4K5B','X-ray 1.85 A','entity2 C/D:184 aa；核心缺Met1，P117V/Q133R；SIFTS Q1RMX3/Bos taurus','entity1 A/B:DARPin；描述名称易误导','受体可结合构象参考','不能按description选链；同样不是精确assay或BH3肽复合物','https://www.rcsb.org/structure/4K5B'),
    ]
    writecsv('structure_source_input_map.csv', [dict(zip(
        ('pdb_id','method','receptor_correspondence','partner_correspondence','allowed_use','exclusion_or_caveat','primary_url'),r)) for r in structures])
    hypotheses = [
        ('H1','single-sequence缺受体进化信息使折叠/对接不充分','六A3M无同源行；受体均值39.92-41.32且无残基>=70；历史成功MDM2使用MSA',
         '跨体系比较非受控；没有Bcl-w同输入MSA对照，不能证明因果','优先解释候选，未证实','未来若有合规共同受体MSA，固定其他项仅改MSA；当前无可冻结MSA，HOLD'),
        ('H2','加速内核/数值实现路径使Bcl-w低置信度','当前启用fast kernels/BF16/fuse；日志有计时与int64截断警告；没有当次受控阳性',
         '两警告本身不证明数值错；历史成功配置同版本且同fast/BF16/fuse；整套环境坏并无证据','次级、可低成本区分','已准备仅亲本双seed的--no-use-fast-kernels单选项诊断；本轮未执行'),
        ('H3','受体构建/天然变体与实验对象或结构参考不一致','assay完整克隆缺失；旧结构多R133，当前Q133；UniProt确认Q133R天然变体',
         '区间1-164与论文方法一致；Q133R影响未知；不能由旧结构倒推assay一定R133','真实输入不确定性，不是已证实输入错误','恢复公开克隆序列/版本后再定；不擅自换变体或延长尾部'),
        ('H4','荧光团/酰胺/标签省略导致化学与结构差异','原文直接Kd对象带染料/C-amide/c-myc，模型无这些组分',
         '缺带/不带修饰同条件结构对照；AFM本来只编码标准核心；不能证明省略造成全受体低置信','明确适用范围限制，因果未知','保留代理定义；无精确接头不凭经验补齐，不直接建化学物救分'),
        ('H5','链顺序、长度、残基索引、FASTA或A3M损坏','专门核查输入、A3M、187x187置信度、PDB A164/B23',
         '六输入/十二PDB一致；PAE导出与分数矩阵一致；265原文件无偏差','本次检查范围内不支持','无需重跑映射或更换输入；保留核验表'),
        ('H6','放松或scorer导致初始置信度失败','原始结构已有低pLDDT与高PAE；放松前存在严重近接；CA位移0.185-0.555 A',
         '放松仅消除几何近接，不重新预测置信度；GraphPep/FoldX是下游','不能解释上游低置信度','不重复放松/评分，不把放松成功当结构正确'),
        ('H7','冻结PAE最大值太严格是唯一失败原因','pilot使用全双向跨链项最大值；老协议写中位数，命名为interface易混',
         '六项ipTM与肽pLDDT也各自失败；跨链中位数21.28-23.22 A仍>10','不能作为唯一原因或事后放宽理由','旧阈值不动；新方案沿用，所有替代PAE统计只描述'),
        ('H8','文献肽不结合所以预测低置信','本模型低置信且seed姿态不稳定',
         '存在来源直接Kd；置信度不是实验结合试验；没有新的生物学测量','不允许该结论','保留文献标签语义，不据此宣布不结合'),
    ]
    writecsv('root_cause_hypotheses_evidence.csv', [dict(zip(
        ('id','hypothesis','support','counterevidence_or_missing_control','assessment','bounded_next_step'),r)) for r in hypotheses])
    prepared = OUT / 'prepared_inputs'
    prepared.mkdir(exist_ok=True)
    lines = (PILOT / 'complexes.fasta').read_text().splitlines()
    parent_fasta = prepared / 'parent_only.fasta'
    parent_fasta.write_text('\n'.join(lines[:2])+'\n',encoding='utf-8')
    frozen = json.loads((PILOT / 'run_config_frozen_v1.json').read_text())
    run_root = ROOT / 'output/bclw_parent_kernel_ablation_v1'
    assert not run_root.exists(), 'Review existing run rather than overwrite'
    commands = []
    for seed in (20261008,20261009):
        event = next(e for e in json.loads((PILOT / 'run_execution_record_v1.json').read_text())['events']
                     if e['label']==f'afm_seed_{seed}')
        argv = list(event['command'])
        argv[1] = str(parent_fasta)
        argv[2] = str(run_root/f'seed_{seed}/structures')
        argv[argv.index('--use-fast-kernels')] = '--no-use-fast-kernels'
        commands.append({'seed':seed,'environment':{'XLA_PYTHON_CLIENT_PREALLOCATE':'false'},
                         'argv':argv,'relax_argv':[frozen['prediction']['python'],str(ROOT/'pipeline/relax_feedback_panel.py'),
                         str(run_root/f'seed_{seed}/structures'),str(run_root/f'seed_{seed}/relaxed')]})
    plan = {
        'protocol_id':'BCLW_PARENT_KERNEL_ABLATION_V1',
        'status':'PREPARED_NOT_AUTHORIZED_TO_RUN_IN_THIS_TURN_NOT_STARTED',
        'question':'With identical core sequences, seeds, model, no MSA/templates and recycles, is disabling the fast-kernels path alone sufficient to recover confident parent predictions?',
        'purpose':'Implementation-path diagnostic only; not a candidate-reward comparison or a source of affinity labels.',
        'why_small_compute_may_be_worthwhile':'A two-parent branch can separate one implementation factor without MSA service, new data selection, chemical modeling or six mutant re-runs. Historical fast-kernel success lowers prior suspicion; this is a cheap diagnostic, not a presumed fix.',
        'authorization':'This turn prepares only. A later explicit user instruction to start this prepared protocol authorizes its execution; no repeat confirmation required within these bounds.',
        'input_fasta':str(parent_fasta),'input_fasta_sha256':sha(parent_fasta),
        'input_chemistry':'Unchanged Q92843 1-164 + standard Bim_A2eT 23mer; core proxy, Q133 and all omitted chemistry retained.',
        'baseline':{'directory':str(PILOT),'parent_seeds':[20261008,20261009],
                    'frozen_config_sha256':sha(PILOT/'run_config_frozen_v1.json'),
                    'reuse':'Compare with existing same-seed parent outputs; do not rerun baseline or rescore it.'},
        'single_changed_prediction_option':'--use-fast-kernels -> --no-use-fast-kernels',
        'constant_prediction_settings':{key:frozen['prediction'][key] for key in (
            'python','cli','versions','model_type','model_number','model_parameter_file','model_parameter_sha256',
            'templates','recycles','recycle_early_stop_tolerance','max_seq','max_extra_seq','num_ensemble',
            'pair_mode','num_relax_in_colabfold','seeds','msa_policy')},
        'constant_relax_settings':frozen['prediction']['commands']['relax_settings'],
        'commands':commands,'planned_output_root':str(run_root),
        'budget':{'new_prediction_jobs_max':2,'raw_final_structures_max':2,'restrained_relaxed_structures_max':2,
                  'gpu_tagged_process_wall_seconds_max':1800,'cpu_relax_wall_seconds_max':900,'disk_bytes_max':536870912,
                  'run_serially':True,'budget_meaning':'Process wall-time cap, not integrated GPU utilization.',
                  'no_reward_calls':True,'no_msa_or_template_service':True,'no_weight_download_or_reinstall':True},
        'preflight_before_future_start':[
            'Re-read applicable AGENTS and latest user steering; confirm original integrity and input/config SHA256.',
            'Check executable, installed batch source SHA256 a65a1fb4301ae82de8b307a8d698f327a79b8fa66c8f053f4ddeff52c87020c4 and model weights against original hashes; no version switch.',
            'Record CUDA/JAX/runtime versions and GPU free memory; GPU total about 8GB is not a reservation; do not kill unrelated jobs.',
            'Confirm output root absent; if occupied prepare and record unique new root before any prediction.',
            'Use a budget watchdog for process wall/disk caps; stop on OOM, NaN, error or changed inputs; retain all failures, do not change seed/model/recycles to rescue.',
            'Retain baseline effective MSA1/1 and all derived config fields; report differences derived from the kernel option rather than concealing them.'
        ],
        'quality_gate_unchanged':frozen['proposed_readouts_pre_score']['quality_gate'],
        'other_required_gates':['Exact A164/B23 sequences and contiguous numbering','At least one unique heavy-atom-defined residue contact <=6 A','After identical restrained Amber relaxation, minimum interchain heavy-atom distance >=1.5 A'],
        'additional_diagnostic_gate':{'receptor_mean_plddt_min':70,'note':'Added engineering readiness check; never retroactively changes original pilot pass count.'},
        'frozen_branch_decision':{
            '2_of_2_all_gates_pass':'Evidence that the changed implementation path is sufficient for this parent under these seeds; consider a separately prepared continuation, with mutant scoring still unauthorized.',
            '1_of_2_all_gates_pass':'Seed-sensitive/partial result; HOLD. No reroll or claim of a rescued reproducible system.',
            '0_of_2_all_gates_pass':'Disabling fast kernels alone does not recover the parent. HOLD this branch; does not prove all runtime factors innocent or MSA/chemistry causal.',
            'confidence_improves_without_gate_pass':'Report only descriptive change, classify as gate failure; no new success threshold.',
            'any_execution_failure':'Count failed job and stop; preserve logs; no automatic fallback or new parameters.'},
        'excluded':['GraphPep','FoldX','ESM2','PepMLM','mutant predictions','full panel','MD/P27A minimization','training','chemical reconstruction','wet lab','paid services'],
        'source_and_direction_boundaries':{'full_benchmark':'HOLD; insufficient source labels/backgrounds and exact assay input',
             'given_pose_scoring':'HOLD; any future template/experimental-pose branch changes the scientific question and requires separate pose/mapping freeze.',
             'future_ties_and_denominators':'Do not change v2 0.05 GraphPep/FoldX tie bands; all formally eligible labels retain failures/ties in denominator; not calculated here.'}
    }
    writej('next_compute_protocol_prepared_v1.json',plan)
    writej('decision.json',{'audit_status':'COMPLETE_WITH_CAUSAL_UNCERTAINTY','new_prediction_jobs_started':0,'new_reward_calls':0,
        'historical_pilot':'COMPLETED; 0/6 frozen quality gate; not zero direction accuracy',
        'go':'GO_TO_REVIEW_PREPARED_TWO_PARENT_KERNEL_DIAGNOSTIC_ONLY',
        'execution':'HOLD_THIS_TURN; user authorized preparation and read-only audit only',
        'msa_branch':'HOLD_NO_AUDITED_RECEPTOR_HOMOLOG_MSA_INPUT; do not fabricate homolog rows from six pilot A3Ms',
        'experimental_pose_branch':'HOLD_NO_EXACT_ASSAY_COMPLEX; related structures must not be relabeled as Bcl-w:Bim_A2eT',
        'full_benchmark_and_design':'HOLD; data/structure/reward gates not passed',
        'root_cause':'Not established; limited evidence disfavors mapping corruption and downstream scoring as causes of initial low confidence.',
        'unresolved':['MSA effect','target-specific numerical/runtime effects','Q133/R133 assay clone','full tag/linker','chemical/aggregation effects','full independent bidirectional benchmark']})
    writej('audit_scope_and_rules.json',{'date':'2026-10-20 Asia/Shanghai','workspace':str(ROOT),
        'handoff_read':'PROJECT_HANDOFF_2026-10-20.md, full',
        'rules':{'E:/AGENTS.md':'Read; FastAPI-specific. No FastAPI application or API modifications in this research audit.',
                 'E:/Users/AGENTS.md':'absent','workspace/AGENTS.md':'absent','workspace_descendants':'rg --files --hidden -g AGENTS.md found no file',
                 'wsl_package_ancestors':'Checked /, /home, /home/g and colabfold batch.py ancestor chain; no AGENTS found',
                 'wsl_mount_parent':'/mnt/AGENTS.md absent; /mnt/e/AGENTS.md maps to already-read E:/AGENTS.md'},
        'skills':['systematic-debugging','verification-before-completion','python-dev','pdf read-only source workflow'],
        'mutation_boundary':'Only new files under this audit directory; no original outputs, frozen files, package installations or long-term memories changed.',
        'analysis_boundary':'Small local parsing, hashing and confidence/geometric audit of existing structures are not new model/reward inference.',
        'database_scope':'5 Q92843 accession-linked RCSB entries + 3FDL + 4K5A/4K5B primary-paper references; 8 unique entries. Not an exhaustive proof of nonexistence.',
        'retrieval_failure':'EuropePMC PMC4276867 XML HTTP500 retained as failed; primary PMC HTML and RCSB data obtained. No Dutta2013 supplement redownload.',
        'memory_use':'Old memory used only for boundary/path orientation; current handoff and raw evidence supersede old HOLD and not-started fields.'})
    print('Prepared tables, parent FASTA, decision and protocol. No inference or scoring started.')


if __name__ == '__main__':
    main()
