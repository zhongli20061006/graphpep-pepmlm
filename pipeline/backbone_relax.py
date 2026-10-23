"""Common hydration and isolated OpenMM minimization for fixed/flexible arms.

No dynamics, imposed helicity, or use of potential energy as binding energy.
Imports OpenMM lazily so preparation and raw scoring audits need no OpenMM.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import random
import re
import numpy as np
from pipeline import backbone_geometry as geo
from pipeline.foldx_matched_diagnostic import require, sha, save_json, load_json, validate_pdb, check_hashes

FORCEFIELDS=('amber14/protein.ff14SB.xml','implicit/obc2.xml')


def frozen_indices(atoms: list, condition: str) -> list[int]:
    require(condition in ('fixed','flexible'),'Unknown relaxation condition')
    return [i for i,(chain,_,name) in enumerate(atoms) if name in geo.BACKBONE and (chain=='A' or (chain=='B' and condition=='fixed'))]


def tether_indices(atoms: list, condition: str) -> list[int]:
    require(condition in ('fixed','flexible'),'Unknown relaxation condition')
    return [i for i,(chain,_,name) in enumerate(atoms) if chain=='B' and name in geo.BACKBONE] if condition=='flexible' else []


def tether_k(kcal_mol_A2: float) -> float:
    return kcal_mol_A2*4.184*100


def strip_hydrogens(text: str) -> str:
    lines=[]
    for line in text.splitlines():
        if line.startswith(('ATOM  ','HETATM')):
            name=line[12:16].strip()
            element=(line[76:78].strip() if len(line)>=78 else '') or re.sub(r'^[0-9]+','',name)[0]
            if element.upper() in ('H','D'): continue
        lines.append(line)
    return '\n'.join(lines)+'\n'


def atom_mapping(topology) -> list[dict]:
    return [dict(index=a.index,chain=a.residue.chain.id,residue=int(a.residue.id),residue_name=a.residue.name,atom=a.name,element=a.element.symbol) for a in topology.atoms()]


def mapping_keys(mapping: list[dict]) -> list[tuple]:
    return [(a['chain'],a['residue'],a['atom']) for a in mapping]


def forcefield():
    from openmm import app
    return app.ForceField(*FORCEFIELDS)


def hydrate(source: Path, common: Path, source_uid: str, sequence: str) -> dict:
    """Exactly one seeded H preparation per immutable source identity."""
    from openmm import app, unit, Platform
    import io
    source,common=Path(source),Path(common)
    source_hash=sha(source)
    if (common/'hydration.json').exists():
        metadata=load_json(common/'hydration.json')
        require(metadata['source_uid']==source_uid and metadata['source_sha256']==source_hash,'Common hydration source changed')
        check_hashes(metadata['artifact_hashes'],common)
        return metadata
    require(not common.exists(),'Untracked/incomplete common hydration; preserve for inspection')
    common.mkdir(parents=True)
    stripped=strip_hydrogens(source.read_text(encoding='ascii'))
    validate_pdb(stripped,sequence)
    original=geo.coordinates(stripped)
    (common/'heavy_start.pdb').write_text(stripped,encoding='ascii')
    pdb=app.PDBFile(io.StringIO(stripped)); modeller=app.Modeller(pdb.topology,pdb.positions)
    digest=hashlib.sha256(f'20261007:{source_uid}:{source_hash}'.encode()).digest()
    seed=int.from_bytes(digest[:8],'big')%(2**32)
    random.seed(seed)
    platform=Platform.getPlatformByName('CPU')
    properties={'Threads':'1','DeterministicForces':'true'}
    for key,value in properties.items(): platform.setPropertyDefaultValue(key,value)
    modeller.addHydrogens(forcefield(),pH=7.4,platform=platform)
    mapping=atom_mapping(modeller.topology)
    positions=np.asarray(modeller.positions.value_in_unit(unit.angstrom),dtype=float)
    heavy={k:tuple(positions[i]) for i,k in enumerate(mapping_keys(mapping)) if mapping[i]['element'] not in ('H','D')}
    delta=geo.displacement(original,heavy)
    require(delta['max_A']<=1e-9,'Hydration moved heavy atoms')
    require(len(mapping)==len(set(mapping_keys(mapping))),'Ambiguous hydrated atom identities')
    with (common/'hydrated.pdb').open('w') as stream:
        app.PDBFile.writeFile(modeller.topology,modeller.positions,stream,keepIds=True)
    save_json(common/'positions.json',dict(units='angstrom',atoms=mapping,positions=positions.tolist()))
    metadata=dict(source_uid=source_uid,source_sha256=source_hash,hydrogen_seed=seed,seed_derivation='sha256(20261007:source_uid:source_sha256) first 8 bytes big-endian modulo 2**32',heavy_max_displacement_A=delta['max_A'],pH=7.4,platform='CPU',platform_properties=properties,artifact_hashes={p.name:sha(p) for p in common.iterdir() if p.is_file()})
    save_json(common/'hydration.json',metadata)
    return metadata


def create_system(topology, start_nm: np.ndarray, mapping: list[dict], condition: str):
    from openmm import app, unit, CustomExternalForce
    system=forcefield().createSystem(topology,nonbondedMethod=app.NoCutoff,constraints=None,rigidWater=False,removeCMMotion=False,soluteDielectric=1,solventDielectric=78.5,implicitSolventKappa=1.274/unit.nanometer)
    atoms=mapping_keys(mapping); frozen=frozen_indices(atoms,condition)
    for i in frozen: system.setParticleMass(i,0*unit.dalton)
    for force in system.getForces(): force.setForceGroup(0)
    tether=CustomExternalForce('0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    tether.addGlobalParameter('k',41.84)
    for name in ('x0','y0','z0'): tether.addPerParticleParameter(name)
    for i in tether_indices(atoms,condition): tether.addParticle(i,start_nm[i].tolist())
    tether.setForceGroup(1); system.addForce(tether)
    return system,frozen


def force_metrics(forces, frozen: list[int]) -> dict:
    arr=np.asarray(forces,dtype=float); free=np.ones(len(arr),dtype=bool); free[frozen]=False
    return dict(all_force_rms_kJ_mol_nm=float(np.sqrt(np.mean(np.sum(arr**2,axis=1)))),free_force_rms_kJ_mol_nm=float(np.sqrt(np.mean(np.sum(arr[free]**2,axis=1)))),free_force_component_rms_kJ_mol_nm=float(np.sqrt(np.mean(arr[free]**2))),masked_gradient_component_rms_kJ_mol_nm=float(np.sqrt(np.sum(arr[free]**2)/(3*len(arr)))))


def relax(common: Path, output: Path, condition: str, sequence: str, max_iterations=2000, platform_name='CPU') -> dict:
    from openmm import app, unit, Context, Platform, VerletIntegrator, LocalEnergyMinimizer, MinimizationReporter
    common,output=Path(common),Path(output)
    metadata=load_json(common/'hydration.json'); check_hashes(metadata['artifact_hashes'],common)
    data=load_json(common/'positions.json'); mapping=data['atoms']; start_A=np.asarray(data['positions'],dtype=float)
    pdb=app.PDBFile(str(common/'hydrated.pdb'))
    require(atom_mapping(pdb.topology)==mapping,'Hydrated PDB topology disagrees with full precision identity mapping')
    start_nm=start_A/10
    system,frozen=create_system(pdb.topology,start_nm,mapping,condition)
    integrator=VerletIntegrator(.001*unit.picosecond)
    properties={'Threads':'1','DeterministicForces':'true'} if platform_name=='CPU' else {}
    context=Context(system,integrator,Platform.getPlatformByName(platform_name),properties)
    context.setPositions(start_nm*unit.nanometer)
    def snapshot():
        state=context.getState(getEnergy=True,getForces=True)
        metrics=force_metrics(state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/unit.nanometer),frozen)
        metrics['physical_potential_kJ_mol']=context.getState(getEnergy=True,groups={0}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        metrics['tether_potential_kJ_mol']=context.getState(getEnergy=True,groups={1}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        require(all(np.isfinite(v) for v in metrics.values()),'Nonfinite minimization energy/force')
        return metrics
    before=snapshot()
    class Reporter(MinimizationReporter):
        def __init__(self): super().__init__(); self.calls=[]
        def report(self,iteration,x,grad,args):
            self.calls.append(dict(iteration=int(iteration),system_energy_kJ_mol=float(args['system energy']),restraint_energy_kJ_mol=float(args['restraint energy']),max_constraint_error=float(args['max constraint error'])))
            return False
    reporter=Reporter()
    LocalEnergyMinimizer.minimize(context,10*unit.kilojoule_per_mole/unit.nanometer,int(max_iterations),reporter)
    after=snapshot()
    final_A=np.asarray(context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom),dtype=float)
    require(np.isfinite(final_A).all(),'Nonfinite minimized positions')
    frozen_max=float(np.linalg.norm(final_A[frozen]-start_A[frozen],axis=1).max())
    require(frozen_max<=1e-6,'Frozen atom moved beyond tolerance before PDB rounding')
    output.mkdir(parents=True,exist_ok=True)
    save_json(output/'positions.json',dict(units='angstrom',atoms=mapping,positions=final_A.tolist(),start_positions_sha256=sha(common/'positions.json')))
    with (output/'relaxed_with_H.pdb').open('w') as stream:
        app.PDBFile.writeFile(pdb.topology,final_A/10*unit.nanometer,stream,keepIds=True)
    final_heavy={k:tuple(final_A[i]) for i,k in enumerate(mapping_keys(mapping)) if mapping[i]['element'] not in ('H','D')}
    start_heavy={k:tuple(start_A[i]) for i,k in enumerate(mapping_keys(mapping)) if mapping[i]['element'] not in ('H','D')}
    lines=[]
    for line in (common/'heavy_start.pdb').read_text().splitlines():
        if line.startswith('ATOM  '):
            key=(line[21],int(line[22:26]),line[12:16].strip()); x,y,z=final_heavy[key]
            require(max(abs(x),abs(y),abs(z))<999,'Coordinate cannot be represented safely in PDB')
            line=line[:30]+f'{x:8.3f}{y:8.3f}{z:8.3f}'+line[54:]
        lines.append(line)
    heavy='\n'.join(lines)+'\n'; validate_pdb(heavy,sequence)
    require(geo.coordinates(heavy).keys()==start_heavy.keys(),'Relaxed heavy identity changed')
    (output/'relaxed.pdb').write_text(heavy,encoding='ascii')
    result=dict(condition=condition,hydration_sha256=sha(common/'hydration.json'),start_positions_sha256=sha(common/'positions.json'),platform=platform_name,platform_properties=properties,iterations=len(reporter.calls),last_reported_iteration=reporter.calls[-1]['iteration'] if reporter.calls else None,iteration_cap=int(max_iterations),iteration_cap_reached=len(reporter.calls)>=max_iterations,converged=after['free_force_rms_kJ_mol_nm']<=10,convergence_definition='free atom vector-force RMS <= 10 kJ/mol/nm; iteration cap is separately recorded',frozen_atom_indices=frozen,frozen_max_A=frozen_max,before=before,after=after,reporter=reporter.calls,**geo.movement(start_heavy,final_heavy),**geo.geometry(final_heavy,sequence))
    result['iteration_measurement']='MinimizationReporter call count; not proof of optimizer stopping criterion'
    result['torsions_before']=geo.torsions(start_heavy,sequence)
    result['torsions_after']=geo.torsions(final_heavy,sequence)
    result['hbond_proxies_before']=geo.hbond_proxies(start_heavy,sequence)
    result['hbond_proxies_after']=geo.hbond_proxies(final_heavy,sequence)
    result['residue_displacements']=[]
    for number in range(1,len(sequence)+1):
        row=dict(residue=number,amino_acid=sequence[number-1])
        for bb,prefix in [(True,'backbone'),(False,'all_heavy')]:
            row.update({prefix+'_'+k:v for k,v in geo.displacement(geo.select(start_heavy,'B',bb,number,number),geo.select(final_heavy,'B',bb,number,number)).items()})
        result['residue_displacements'].append(row)
    save_json(output/'relaxation.json',result)
    del context,integrator
    return result
