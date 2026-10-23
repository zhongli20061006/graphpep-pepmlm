"""Behavioral tests for the frozen paired backbone sensitivity diagnostic."""
from __future__ import annotations
import importlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class BackboneDiagnosticTests(unittest.TestCase):
    def module(self, name):
        self.assertIsNotNone(importlib.util.find_spec('pipeline.' + name), 'Required diagnostic feature is absent: ' + name)
        return importlib.import_module('pipeline.' + name)

    def test_mass_selection_and_tether_units(self):
        r = self.module('backbone_relax')
        atoms = [('A', 1, 'N'), ('A', 1, 'CB'), ('B', 1, 'CA'), ('B', 1, 'H'), ('B', 12, 'OXT')]
        self.assertEqual(r.frozen_indices(atoms, 'fixed'), [0, 2, 4])
        self.assertEqual(r.frozen_indices(atoms, 'flexible'), [0])
        self.assertEqual(r.tether_indices(atoms, 'flexible'), [2, 4])
        self.assertEqual(r.tether_indices(atoms, 'fixed'), [])
        self.assertAlmostEqual(r.tether_k(0.1), 41.84)
        with self.assertRaises(ValueError):
            r.frozen_indices(atoms, 'unfrozen')

    def test_direct_displacement_and_kabsch_receptor_alignment(self):
        g = self.module('backbone_geometry')
        a = {('A', 1, 'CA'): (0., 0., 0.), ('A', 2, 'CA'): (1., 0., 0.), ('A', 3, 'CA'): (0., 1., 0.), ('B', 1, 'CA'): (0., 0., 1.)}
        b = {k: (-v[1] + 3, v[0] + 4, v[2] + 5) for k, v in a.items()}
        self.assertGreater(g.displacement(a, b)['rmsd_A'], 1)
        aligned = g.align_on_receptor(a, b)
        self.assertLess(g.displacement(a, aligned)['max_A'], 1e-10)
        bad = dict(a); bad.pop(('B', 1, 'CA'))
        with self.assertRaises(ValueError):
            g.displacement(a, bad)

    def test_torsion_and_terminal_missingness_and_proline_donor(self):
        g = self.module('backbone_geometry')
        self.assertAlmostEqual(abs(g.dihedral((1,0,0), (0,0,0), (0,1,0), (0,1,1))), 90)
        coords = {}
        for i in range(1, 13):
            for atom, xyz in [('N',(i,0,0)), ('CA',(i,.5,.5)), ('C',(i,1,0)), ('O',(i+.1,1,0))]:
                coords[('B',i,atom)] = xyz
        rows = g.torsions(coords, 'ETFSDLWKLLPE')
        self.assertIsNone(rows[0]['phi_deg']); self.assertIsNone(rows[-1]['psi_deg'])
        self.assertEqual(len([x for x in rows if x['interior_2_11']]), 10)
        proxies = g.hbond_proxies(coords, 'ETFSDLWKLLPE')
        self.assertFalse(next(x for x in proxies if x['donor_residue'] == 11)['eligible_donor'])
        self.assertFalse(next(x for x in proxies if x['donor_residue'] == 11)['proxy_present'])

    def test_paired_contrast_is_median_of_corresponding_differences(self):
        a = self.module('backbone_analysis')
        rows = []
        for rep, fixed, flexible in [(1,0,100),(2,100,101),(3,101,0),(4,102,103),(5,103,104)]:
            for condition, gain in [('fixed',fixed),('flexible',flexible)]:
                rows.append(dict(pair_uid=f'x_r{rep}', peptide_id='P53_P27A', seed=20261011, model_number=1, repeat=rep, condition=condition, status='ok', gain_kcal_mol=gain))
        contrasts = a.paired_contrasts(rows)
        self.assertEqual(len(contrasts), 5)
        summary = a.contrast_summary(contrasts)
        self.assertEqual(summary[0]['median_flexible_minus_fixed_gain_kcal_mol'], 1)
        self.assertNotEqual(summary[0]['median_flexible_minus_fixed_gain_kcal_mol'], 101-101)
        rows.pop()
        missing = a.paired_contrasts(rows)
        self.assertEqual(next(x for x in missing if x['repeat'] == 5)['status'], 'missing')
        self.assertEqual(a.contrast_summary(missing)[0]['status'], 'incomplete')

    def test_technical_gate_score_agnostic_and_partial_never_passes(self):
        a = self.module('backbone_analysis')
        rows = [dict(job_id=str(i), status='ok', condition='flexible' if i%2 else 'fixed', frozen_max_A=0., B_backbone_max_A=.1 if i%2 else 0., converged=False, energy_kcal_mol=1000*(-1)**i) for i in range(48)]
        self.assertTrue(a.technical_gate(rows, 48)['technical_pass'])
        self.assertFalse(a.technical_gate(rows[:-1], 48)['technical_pass'])
        rows[0]['status']='failed'
        self.assertFalse(a.technical_gate(rows,48)['technical_pass'])

    def test_strict_source_identity_rejects_wrong_mutation_chain_or_oxt(self):
        d = self.module('backbone_diagnostic')
        source = ROOT/'output/foldx_matched_diagnostic_20261006/jobs/P53_P27A_seed20261011_model1/input_1_0.pdb'
        text = source.read_text()
        d.validate_source(text, 'ETFSDLWKLLAE')
        with self.assertRaises(ValueError): d.validate_source(text, 'ETFSDLWKLLPE')
        with self.assertRaises(ValueError): d.validate_source('\n'.join(x for x in text.splitlines() if x[12:16].strip()!='OXT'), 'ETFSDLWKLLAE')
        bad = '\n'.join(x[:21]+'C'+x[22:] if x.startswith('ATOM') and x[21]=='B' else x for x in text.splitlines())
        with self.assertRaises(ValueError): d.validate_source(bad, 'ETFSDLWKLLAE')

    def test_summary_required_and_intraclashes_valid(self):
        d = self.module('backbone_diagnostic')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            raw='Pdb\tGroup1\tGroup2\tInteraction Energy\tIntraclashesGroup1\tIntraclashesGroup2\nx.pdb\tA\tB\t-1\t2\t3\n'
            (p/'Interaction_x.fxout').write_text(raw)
            with self.assertRaises(ValueError): d.parse_score(p, 'x.pdb')
            (p/'Summary_x.fxout').write_text(raw)
            self.assertEqual(d.parse_score(p,'x.pdb')['energy_kcal_mol'], -1)
            (p/'Summary_x.fxout').write_text(raw.replace('-1','-2'))
            with self.assertRaises(ValueError): d.parse_score(p,'x.pdb')

    def test_preserve_interrupted_attempt_and_changed_hash_is_permanent_failure(self):
        d = self.module('backbone_diagnostic')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); (p/'stdout.txt').write_text('interrupted evidence')
            archived=d.archive_attempt(p)
            self.assertEqual((archived/'stdout.txt').read_text(),'interrupted evidence')
            (p/'raw.txt').write_text('raw')
            result=dict(status='failed',job_id='x',signature='sig',failure_reason='process failed',provenance_valid=True,artifact_hashes={'raw.txt':d.sha(p/'raw.txt')})
            d.save_json(p/'result.json',result)
            expected=d.sha(p/'result.json')
            self.assertEqual(d.validate_cached(p,'x','sig',expected)['status'],'failed')
            (p/'raw.txt').write_text('changed')
            with self.assertRaises(ValueError): d.validate_cached(p,'x','sig',expected)
            d.invalidate_job(p,'x','sig','hash changed')
            self.assertFalse(d.load_json(p/'result.json')['provenance_valid'])
            with self.assertRaises(ValueError): d.require_promotable(p)

    @unittest.skipUnless(importlib.util.find_spec('openmm'), 'OpenMM integration needs frozen WSL environment')
    def test_real_openmm_common_hydration_frozen_and_movable_atoms(self):
        r=self.module('backbone_relax')
        source=ROOT/'output/foldx_matched_diagnostic_20261006/jobs/P53_P27A_seed20261011_model1/input_1_0.pdb'
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            r.hydrate(source,p/'common','integration','ETFSDLWKLLAE')
            d=self.module('backbone_diagnostic')
            self.assertTrue(hasattr(d,'recover_untracked_hydration'),'Completed worker hydration needs strict recovery before parent-state persistence')
            row=dict(source_uid='integration',source_sha256=d.sha(source),sequence='ETFSDLWKLLAE')
            recovered=d.recover_untracked_hydration(source,p/'common',row)
            self.assertEqual(recovered,d.sha(p/'common'/'hydration.json'))
            metadata=d.load_json(p/'common'/'hydration.json')
            metadata['source_uid']='wrong'
            d.save_json(p/'common'/'hydration.json',metadata)
            with self.assertRaises(ValueError): d.recover_untracked_hydration(source,p/'common',row)
            metadata['source_uid']='integration'; d.save_json(p/'common'/'hydration.json',metadata)
            first=r.relax(p/'common',p/'fixed','fixed','ETFSDLWKLLAE',max_iterations=3)
            second=r.relax(p/'common',p/'flexible','flexible','ETFSDLWKLLAE',max_iterations=3)
            self.assertLessEqual(first['frozen_max_A'],1e-6)
            self.assertLessEqual(second['frozen_max_A'],1e-6)
            self.assertGreater(second['B_backbone_max_A'],.001)
            self.assertEqual(first['hydration_sha256'],second['hydration_sha256'])
            self.assertFalse(first['converged'])
            self.assertIn('torsions_before', first)
            self.assertIn('torsions_after', second)
            self.assertEqual(len(second['residue_displacements']),12)
            self.assertEqual(r.load_json(p/'common'/'hydration.json')['platform_properties'], {'Threads':'1','DeterministicForces':'true'})
            print('Temporary OpenMM integration:', json.dumps({condition:{key:result[key] for key in ('iterations','frozen_max_A','B_backbone_max_A','converged')} for condition,result in [('fixed',first),('flexible',second)]},sort_keys=True))
            (p/'common'/'hydrated.pdb').write_text('changed artifact')
            with self.assertRaises(ValueError): d.recover_untracked_hydration(source,p/'common',row)

    def test_prepare_is_immutable_and_complete_without_scoring(self):
        d=self.module('backbone_diagnostic')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'protocol.json').write_bytes((ROOT/'output/backbone_diagnostic_20261007/protocol.json').read_bytes())
            prepared=d.prepare(p,ROOT)
            self.assertEqual(prepared['source_count'],120)
            self.assertEqual(len(d.read_csv(p/'job_manifest.csv')),240)
            self.assertFalse((p/'jobs').exists())
            self.assertEqual(d.prepare(p,ROOT),prepared)
            baseline=d.analysis.baseline(p,ROOT,d.load_json(p/'protocol.json'))
            self.assertEqual(baseline['independent_structures'],18)
            terms=d.read_csv(p/'analysis'/'baseline_energy_term_contrasts.csv')
            self.assertEqual(len({(x['job_id'],x['repeat']) for x in terms}),60)
            self.assertTrue(all('wt_minus_mutant' in x and 'wt_minus_mutant_kcal_mol' not in x for x in terms))
            self.assertEqual(next(x['units'] for x in terms if x['term']=='Number of Residues'),'count')
            analyzed=d.analysis.analyze(p,d.read_csv(p/'job_manifest.csv'),{},'full')
            self.assertFalse(analyzed['complete'])
            self.assertEqual(len(d.read_csv(p/'analysis'/'pair_results.csv')),120)
            self.assertEqual(len(d.read_csv(p/'analysis'/'contrasts.csv')),60)
            self.assertEqual(len(d.read_csv(p/'analysis'/'backbone_summary.csv')),24)
            self.assertEqual(len(d.read_csv(p/'analysis'/'seed_summary.csv')),8)
            self.assertTrue((p/'analysis'/'baseline_geometry.csv').exists())
            src=d.read_csv(p/'source_manifest.csv')[0]
            (p/src['input_path']).write_text('modified')
            with self.assertRaises(ValueError): d.prepare(p,ROOT)

    def test_actual_process_failure_preserves_stdout_stderr_and_return_code(self):
        import sys
        d=self.module('backbone_diagnostic')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            actual=d.command_log([sys.executable,'-c','import sys; print("evidence"); print("failure",file=sys.stderr); sys.exit(7)'],p,'fake_worker',30)
            self.assertEqual(actual['returncode'],7)
            self.assertEqual((p/'fake_worker.stdout.txt').read_text().strip(),'evidence')
            self.assertEqual((p/'fake_worker.stderr.txt').read_text().strip(),'failure')
            self.assertGreaterEqual(actual['seconds'],0)

    def test_publish_only_owned_artifacts_and_preserve_unrelated_analysis(self):
        a=self.module('backbone_analysis')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); (p/'analysis').mkdir()
            (p/'analysis'/'baseline.json').write_text('{}')
            (p/'analysis'/'independent_audit.json').write_text('mutable independently owned audit')
            a.publish_hashes(p)
            published=json.loads((p/'analysis'/'published_artifact_hashes.json').read_text())
            self.assertIn('analysis/baseline.json',published)
            self.assertNotIn('analysis/independent_audit.json',published)

    def test_foldx_uses_proven_config_argument_and_frozen_chemistry(self):
        d=self.module('backbone_diagnostic')
        self.assertTrue(hasattr(d,'analyse_command'),'FoldX command construction must be explicit and testable')
        self.assertEqual(d.analyse_command(Path('/tools/foldx')),['/tools/foldx','-f','analyse.cfg'])
        self.assertTrue(hasattr(d,'ANALYSE_CONFIG'))
        self.assertIn('command=AnalyseComplex\n',d.ANALYSE_CONFIG)
        self.assertIn('analyseComplexChains=A,B\n',d.ANALYSE_CONFIG)
        self.assertIn('noNterm=false\nnoCterm=false\n',d.ANALYSE_CONFIG)

    def test_interrupted_hydration_is_preserved_before_same_source_retry(self):
        d=self.module('backbone_diagnostic')
        self.assertTrue(hasattr(d,'archive_incomplete_hydration'),'Incomplete common preparation needs evidence-preserving recovery')
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp); (source/'common').mkdir()
            (source/'common'/'heavy_start.pdb').write_text('partial preparation')
            archived=d.archive_incomplete_hydration(source)
            self.assertEqual((archived/'heavy_start.pdb').read_text(),'partial preparation')
            self.assertFalse((source/'common').exists())
            (source/'common').mkdir(); (source/'common'/'hydration.json').write_text('{}')
            with self.assertRaises(ValueError): d.archive_incomplete_hydration(source)

    def test_foldx_count_terms_keep_values_and_count_units(self):
        a=self.module('backbone_analysis')
        count_names=['Number of Residues','Interface Residues','Interface Residues Clashing','Interface Residues VdW Clashing','Interface Residues BB Clashing']
        values={name:i+1. for i,name in enumerate(count_names)}|{'Interaction Energy':-12.5,'Backbone Hbond':-1.5}
        jobs=[]; results={}
        for role in ('wt_reference','mutant'):
            job=dict(job_id=role,source_uid=role,pair_uid='pair',peptide_id='P53_P27A',seed='20261011',model_number='1',repeat='1',role=role,condition='fixed',experimental_direction='stronger')
            jobs.append(job)
            results[role]=dict(status='ok',provenance_valid=True,raw=dict(energy_kcal_mol=-12.5,tables={'Interaction_relaxed.fxout':values}),relaxation=dict(before={},after={}))
        _,_,terms=a.assemble_pairs(jobs,results)
        for row in terms:
            self.assertEqual(row['value'],values[row['term']])
            self.assertEqual(row['units'],'count' if row['term'] in count_names else 'kcal/mol')


if __name__=='__main__': unittest.main()
