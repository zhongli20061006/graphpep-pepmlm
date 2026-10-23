"""Scientific guards. Arbitrary fixtures live only in TemporaryDirectory."""
from pathlib import Path
import contextlib
import io
import json
import math
import shutil
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
try:
    import foldx_matched_diagnostic as fx
except ModuleNotFoundError:
    fx = None

AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))


def pdb(sequence, receptor, offset=0.):
    """Arbitrary coordinates, never molecular predictions or scored real results."""
    reverse, lines, serial = {v: k for k, v in AA.items()}, [], 1
    for chain, seq in (("A", receptor), ("B", sequence)):
        for n, letter in enumerate(seq, 1):
            for atom, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")) + (("OXT", "O"),) * (n == len(seq)):
                lines.append(f"ATOM  {serial:5d} {atom:^4s} {reverse[letter]:3s} {chain}{n:4d}    "
                             f"{n + offset:8.3f}{serial * .1:8.3f}{serial * -.2:8.3f}{1.:6.2f}{90.:6.2f}          {element:>2s}  ")
                serial += 1
        lines.append("TER")
    return "\n".join(lines + ["END", ""])


def energy(path, name, value=-9.75, chains=("A", "B")):
    path.write_text("FoldX 5 fixture\nPdb\tGroup1\tGroup2\tIntraclashesGroup1\tIntraclashesGroup2\tInteraction Energy\tStability\n"
                    f"./{name}\t{chains[0]}\t{chains[1]}\t1.25\t2.5\t{value}\t-999\n", encoding="utf-8")


def panel(directory):
    root = Path(directory)
    output = root / "matched"
    output.mkdir()
    protocol = json.loads((ROOT / "output/foldx_matched_diagnostic_20261006/protocol.json").read_text(encoding="utf-8"))
    protocol.update(source_root="prior", prior_result_csv="prior/model_results.csv", source_hashes={})
    old = []
    s = protocol["scientific_protocol"]
    for parent, sequence in s["parents"].items():
        for seed in s["seeds"]:
            for model in s["models"]:
                path = root / f"prior/jobs/{parent}_seed{seed}_model{model}/input_Repair.pdb"
                path.parent.mkdir(parents=True)
                path.write_text(pdb(sequence, s["receptor_sequence"]), encoding="ascii")
                protocol["source_hashes"][path.relative_to(root).as_posix()] = fx.sha(path)
                old.append(dict(peptide_id=parent, parent_id=parent, sequence=sequence, seed=seed, model_number=model,
                                status="ok", interaction_energy_kcal_mol=-10., quality_pass="True"))
    for mutation in s["mutations"]:
        for seed in s["seeds"]:
            for model in s["models"]:
                old.append(dict(peptide_id=mutation["peptide_id"], parent_id=mutation["parent_id"], sequence=mutation["sequence"],
                                seed=seed, model_number=model, status="ok", interaction_energy_kcal_mol=-8., quality_pass="True"))
    prior = root / "prior/model_results.csv"
    prior.write_text(fx.csv_text(old), encoding="utf-8")
    protocol["source_hashes"]["prior/model_results.csv"] = fx.sha(prior)
    binary = root / "tools/foldx"
    binary.parent.mkdir()
    binary.write_text("synthetic executable placeholder, never launched", encoding="ascii")
    protocol["source_hashes"]["tools/foldx"] = fx.sha(binary)
    fx.save_json(output / "protocol.json", protocol)
    return root, output, protocol, binary


def build_fixture(job, mutation, receptor):
    files = []
    original = (job / "input.pdb").read_text(encoding="ascii")
    for i in range(5):
        mut, wt = f"input_1_{i}.pdb", f"WT_input_1_{i}.pdb"
        (job / mut).write_text(pdb(mutation["sequence"], receptor), encoding="ascii")
        (job / wt).write_text(original, encoding="ascii")
        files.extend((mut, wt))
    (job / "PdbList_input.fxout").write_text("\n".join(files) + "\n", encoding="ascii")


def fixture_runner(fail=None):
    """External process stand-in; every output remains in the temporary test root."""
    calls = []
    def runner(command, cwd, stage, timeout):
        calls.append((str(cwd), stage))
        metadata = dict(command=command, cwd=str(cwd), started_utc="SYNTHETIC", seconds=0., returncode=0, error="")
        failure = fail is not None and fail(cwd, stage)
        if failure:
            metadata["returncode"] = 1
        elif stage == "build":
            mutation_text = (cwd / "individual_list.txt").read_text().strip()
            mutation = next(m for m in fx.MUTATIONS if m["individual_list"] == mutation_text)
            build_fixture(cwd, mutation, fx.RECEPTOR)
            for i in range(5):
                (cwd / f"WT_input_1_{i}.pdb").write_text(pdb(fx.PARENTS[mutation["parent_id"]], fx.RECEPTOR, i * .01))
                (cwd / f"input_1_{i}.pdb").write_text(pdb(mutation["sequence"], fx.RECEPTOR, i * .01))
            # Deliberately contradictory total-energy data must never enter binding gains.
            (cwd / "Dif_input.fxout").write_text("Pdb\tTotal Energy\nfixture\t99999\n")
        else:
            name = next(line.split("=", 1)[1] for line in (cwd / "analyse.cfg").read_text().splitlines() if line.startswith("pdb="))
            index = int(name.removesuffix(".pdb").rsplit("_", 1)[1])
            value = -10. - index if cwd.name == "wt_reference" else -12. - 2 * index
            energy(cwd / "Interaction_fixture.fxout", name, value)
            energy(cwd / "Summary_fixture.fxout", name, value)
        fx.save_json(cwd / f"{stage}.command.json", metadata)
        (cwd / f"{stage}.stdout.txt").write_text("SYNTHETIC test only")
        (cwd / f"{stage}.stderr.txt").write_text("")
        return metadata
    return runner, calls


class ScientificGuards(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(fx, "Matched diagnostic implementation does not yet exist")

    def test_energy_named_fields_sign_and_intraclashes(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            path = Path(tmp) / "Interaction_fixture.fxout"
            energy(path, "input.pdb", -9.75, ("B", "A"))
            terms = fx.parse_interaction_terms(path, "input.pdb")
            self.assertEqual(terms["Interaction Energy"], -9.75)
            self.assertEqual(terms["IntraclashesGroup1"], 1.25)
            self.assertEqual(terms["IntraclashesGroup2"], 2.5)
            self.assertEqual(fx.matched_gain(-5., -9.75), 4.75)

    def test_wrong_pdb_chain_and_all_nonfinite_terms_rejected(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            path = Path(tmp) / "Interaction_fixture.fxout"
            for name, value, chains in (("wrong.pdb", -3, ("A", "B")), ("../input.pdb", -3, ("A", "B")), ("input.pdb", -3, ("A", "C")),
                                        ("input.pdb", float("nan"), ("A", "B"))):
                energy(path, name, value, chains)
                with self.assertRaises(ValueError):
                    fx.parse_interaction_terms(path, "input.pdb")
            energy(path, "input.pdb")
            path.write_text(path.read_text().replace("\t1.25\t", "\tinf\t"))
            with self.assertRaises(ValueError):
                fx.parse_interaction_terms(path, "input.pdb")

    def test_duplicate_malformed_and_buildmodel_energy_rejected(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            path = Path(tmp) / "Interaction_fixture.fxout"
            energy(path, "input.pdb")
            original = path.read_text()
            for text in (original + original.splitlines()[-1] + "\n", original.replace("\t-999", ""),
                         original.replace("Stability", "Interaction Energy")):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    fx.parse_interaction_terms(path, "input.pdb")
            path = Path(tmp) / "Dif_input.fxout"
            energy(path, "input.pdb")
            with self.assertRaises(ValueError):
                fx.parse_interaction_terms(path, "input.pdb")

    def test_pdblist_maps_each_reference_and_checks_displacement(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, protocol, _ = panel(tmp)
            fx.prepare(output, root)
            row = fx.read_csv(output / "input_manifest.csv")[0]
            job = output / "jobs" / row["job_id"]
            mutation = protocol["scientific_protocol"]["mutations"][0]
            build_fixture(job, mutation, fx.RECEPTOR)
            pairs = fx.parse_pdblist(job / "PdbList_input.fxout", job, mutation)
            self.assertEqual(len(pairs), 5)
            self.assertEqual([p["repeat"] for p in pairs], [1, 2, 3, 4, 5])
            self.assertEqual(len({p["wt_path"] for p in pairs}), 5)
            self.assertEqual(pairs[0]["pair_backbone_rmsd_A"], 0.)

    def test_missing_duplicate_unlisted_and_wrong_sequence_pairs_rejected(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, protocol, _ = panel(tmp)
            fx.prepare(output, root)
            job = output / "jobs" / fx.read_csv(output / "input_manifest.csv")[0]["job_id"]
            mutation = protocol["scientific_protocol"]["mutations"][0]
            build_fixture(job, mutation, fx.RECEPTOR)
            listing = job / "PdbList_input.fxout"
            original = listing.read_text()
            for text in ("\n".join(original.splitlines()[:-1]) + "\n", original + "input_1_0.pdb\n",
                         original.replace("WT_input_1_4.pdb", "../WT_input_1_4.pdb")):
                listing.write_text(text)
                with self.assertRaises(ValueError):
                    fx.parse_pdblist(listing, job, mutation)
            listing.write_text(original)
            (job / "input_1_0.pdb").write_text((job / "input.pdb").read_text())
            with self.assertRaises(ValueError):
                fx.parse_pdblist(listing, job, mutation)

    def test_missing_oxt_bad_residue_number_and_wrong_chain_rejected(self):
        text = pdb("ETFSDLWKLLPE", fx.RECEPTOR)
        for malformed in ("\n".join(line for line in text.splitlines() if "OXT" not in line),
                          text.replace(" B   1 ", " B   0 "), text.replace(" B", " C")):
            with self.assertRaises(ValueError):
                fx.validate_pdb(malformed, "ETFSDLWKLLPE")

    def test_median_of_five_paired_differences_not_difference_of_medians(self):
        wt, mutant = [0., 100., 101., 102., 103.], [99., 200., 100., 20., 10.]
        gains = [fx.matched_gain(a, b) for a, b in zip(wt, mutant)]
        summary = fx.summarize_backbone(gains)
        self.assertEqual(summary["median_gain_kcal_mol"], 1.)
        self.assertNotEqual(summary["median_gain_kcal_mol"], 101. - 99.)
        self.assertEqual(summary["positive_pairs"], 3)
        self.assertEqual(fx.summarize_backbone(gains[:4])["status"], "incomplete")
        self.assertEqual(fx.summarize_backbone(gains[:4])["median_gain_kcal_mol"], "")

    def test_seed_direction_requires_three_complete_backbones(self):
        for values, expected in (([1, 2, 3], "correct"), ([-1, -2, -3], "incorrect"),
                                 ([1, -1, 3], "inconclusive"), ([1, 0, 3], "inconclusive")):
            summary = fx.summarize_seed(values, "stronger")
            self.assertEqual(summary["diagnostic_direction"], expected)
        self.assertEqual(fx.summarize_seed([1, 2], "stronger")["diagnostic_direction"], "not_evaluable")
        self.assertEqual(fx.summarize_seed([-1, -2, -3], "weaker")["diagnostic_direction"], "correct")

    def test_prepare_24_byte_identical_sources_and_mutation_configs(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, protocol, _ = panel(tmp)
            before = fx.sha(output / "protocol.json")
            fx.prepare(output, root)
            rows = fx.read_csv(output / "input_manifest.csv")
            self.assertEqual(len(rows), 24)
            self.assertEqual([rows[i * 6]["peptide_id"] for i in range(4)], [m["peptide_id"] for m in protocol["scientific_protocol"]["mutations"]])
            for row in rows:
                job = output / "jobs" / row["job_id"]
                self.assertEqual((job / "input.pdb").read_bytes(), (root / row["source_path"]).read_bytes())
                self.assertEqual((job / "individual_list.txt").read_text(), row["individual_list"] + "\n")
                self.assertIn("numberOfRuns=5", (job / "build.cfg").read_text())
            self.assertEqual(before, fx.sha(output / "protocol.json"))
            fx.prepare(output, root)

    def test_source_config_and_frozen_protocol_tampering_rejected(self):
        for target in ("source", "config", "protocol"):
            with self.subTest(target=target), tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
                root, output, _, _ = panel(tmp)
                fx.prepare(output, root)
                row = fx.read_csv(output / "input_manifest.csv")[0]
                path = root / row["source_path"] if target == "source" else output / "protocol.json" if target == "protocol" else output / "jobs" / row["job_id"] / "build.cfg"
                path.write_text(path.read_text() + "\n")
                with self.assertRaises(ValueError):
                    fx.prepare(output, root)

    def test_malformed_mutation_rejected_before_preparation(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, protocol, _ = panel(tmp)
            protocol["scientific_protocol"]["mutations"][0]["individual_list"] = "PB27A;"
            fx.save_json(output / "protocol.json", protocol)
            with self.assertRaises(ValueError):
                fx.prepare(output, root)

    def test_no_results_analysis_is_explicitly_not_run(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, _ = panel(tmp)
            decision = fx.analyze(output, root)
            self.assertEqual(decision["job_counts"], {"ok": 0, "failed": 0, "not_run": 24})
            self.assertEqual(decision["pair_counts"]["not_run"], 120)
            self.assertFalse(decision["complete"])
            self.assertEqual(len(fx.read_csv(output / "analysis/pair_results.csv")), 120)
            self.assertEqual(fx.verify(output, root)["published_tables"], "PASS")

    def test_protocol_aggregation_cannot_disagree_with_computed_rule(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, protocol, _ = panel(tmp)
            protocol["scientific_protocol"]["aggregation"] = "select best pair and best backbone"
            fx.save_json(output / "protocol.json", protocol)
            with self.assertRaises(ValueError):
                fx.prepare(output, root)

    def test_report_preserves_zero_and_warns_about_printed_precision(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, _ = panel(tmp)
            decision = fx.analyze(output, root)
            self.assertIn("printed precision", decision.get("microscopic_gain_caution", ""))
            self.assertIn("舍入", (output / "analysis/report_zh.md").read_text(encoding="utf-8"))
            self.assertEqual(fx.summarize_seed([0., 0., 0.], "weaker")["diagnostic_direction"], "inconclusive")

    def test_frozen_wrong_source_sequence_and_changed_executable_rejected(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, protocol, binary = panel(tmp)
            source_name = next(k for k in protocol["source_hashes"] if "P53_WT" in k)
            source = root / source_name
            original = source.read_bytes()
            source.write_text(pdb("ETFSDLWKLLAE", fx.RECEPTOR))
            protocol["source_hashes"][source_name] = fx.sha(source)
            fx.save_json(output / "protocol.json", protocol)
            with self.assertRaises(ValueError):
                fx.prepare(output, root)
            source.write_bytes(original)
            protocol["source_hashes"][source_name] = fx.sha(source)
            fx.save_json(output / "protocol.json", protocol)
            fx.prepare(output, root)
            binary.write_text("changed synthetic executable")
            with self.assertRaises(ValueError):
                fx.verify(output, root)

    def test_summary_and_interaction_disagreement_rejected(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            directory = Path(tmp)
            energy(directory / "Interaction_fixture.fxout", "input.pdb", -10.)
            energy(directory / "Summary_fixture.fxout", "input.pdb", -11.)
            with self.assertRaises(ValueError):
                fx.raw_energy(directory, "input.pdb")

    def test_missing_or_ambiguous_summary_cannot_be_a_successful_score(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            directory = Path(tmp)
            energy(directory / "Interaction_fixture.fxout", "input.pdb", -10.)
            with self.assertRaises(ValueError):
                fx.raw_energy(directory, "input.pdb")
            energy(directory / "Summary_fixture.fxout", "input.pdb", -10.)
            self.assertEqual(fx.raw_energy(directory, "input.pdb")["energy_kcal_mol"], -10.)
            energy(directory / "Summary_second.fxout", "input.pdb", -10.)
            with self.assertRaises(ValueError):
                fx.raw_energy(directory, "input.pdb")

    def test_actual_per_pair_reference_gain_complete_run_resume_and_raw_audit(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, binary = panel(tmp)
            runner, calls = fixture_runner()
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, workers=2, runner=runner)
            self.assertEqual(len(calls), 264)
            decision = fx.analyze(output, root)
            self.assertTrue(decision["complete"])
            self.assertEqual(decision["pair_counts"]["ok"], 120)
            self.assertEqual(decision["analysecomplex_counts"]["ok"], 240)
            pairs = fx.read_csv(output / "analysis/pair_results.csv")
            self.assertEqual([float(p["gain_kcal_mol"]) for p in pairs[:5]], [2., 3., 4., 5., 6.])
            self.assertEqual(float(pairs[2]["wt_repacking_delta_kcal_mol"]), -2.)
            self.assertEqual(float(pairs[2]["old_independent_gain_kcal_mol"]), -2.)
            self.assertAlmostEqual(float(pairs[2]["wt_source_backbone_rmsd_A"]), .02)
            self.assertEqual(float(pairs[2]["pair_backbone_rmsd_A"]), 0.)
            self.assertEqual(fx.verify(output, root)["published_tables"], "PASS")
            before = len(calls)
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, workers=1, runner=runner)
            self.assertEqual(len(calls), before)
            first = output / "jobs" / pairs[0]["job_id"]
            raw = first / "ac/repeat1/wt_reference/Interaction_fixture.fxout"
            raw.write_text(raw.read_text().replace("-10.0", "-11.0"))
            with self.assertRaises(ValueError):
                fx.verify(output, root)

    def test_one_failed_analysecomplex_cannot_pass_complete_and_preserves_other_pairs(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, binary = panel(tmp)
            first_name = "P53_P27A_seed20261011_model1"
            runner, _ = fixture_runner(lambda cwd, stage: stage == "analyse" and first_name in str(cwd) and "repeat2" in str(cwd) and cwd.name == "mutant")
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, runner=runner)
            decision = fx.analyze(output, root)
            self.assertFalse(decision["complete"])
            self.assertEqual(decision["job_counts"], {"ok": 23, "failed": 1, "not_run": 0})
            self.assertEqual(decision["pair_counts"], {"ok": 119, "failed": 1, "not_run": 0})
            self.assertEqual(decision["analysecomplex_counts"], {"ok": 239, "failed": 1, "not_run": 0})
            backbones = fx.read_csv(output / "analysis/backbone_summary.csv")
            self.assertEqual(backbones[0]["median_gain_kcal_mol"], "")
            seeds = fx.read_csv(output / "analysis/seed_summary.csv")
            self.assertEqual(seeds[0]["diagnostic_direction"], "not_evaluable")
            self.assertEqual(fx.verify(output, root)["published_tables"], "PASS")

    def test_interrupted_outputs_archived_and_untracked_result_rejected(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, binary = panel(tmp)
            fx.prepare(output, root)
            first = output / "jobs" / fx.read_csv(output / "input_manifest.csv")[0]["job_id"]
            (first / "Interrupted_output.fxout").write_text("SYNTHETIC interrupted attempt")
            runner, _ = fixture_runner()
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, runner=runner)
            archived = list((first / "attempts").rglob("Interrupted_output.fxout"))
            self.assertEqual(len(archived), 1)
            self.assertEqual(archived[0].read_text(), "SYNTHETIC interrupted attempt")
            state = fx.load_json(output / "execution_state.json")
            del state["jobs"][first.name]
            fx.save_json(output / "execution_state.json", state)
            with self.assertRaises(ValueError):
                fx.run(output, binary, root, runner=runner)

    def test_published_artifact_tampering_rejected(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, _ = panel(tmp)
            fx.analyze(output, root)
            path = output / "analysis/pair_results.csv"
            path.write_text(path.read_text() + "\n")
            with self.assertRaises(ValueError):
                fx.verify(output, root)

    def test_directory_creation_failure_after_valid_build_keeps_five_auditable_pairs(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, binary = panel(tmp)
            fx.prepare(output, root)
            first_name = fx.read_csv(output / "input_manifest.csv")[0]["job_id"]
            runner, calls = fixture_runner()
            original_mkdir = Path.mkdir
            failed = False
            def fail_one_directory(path, *args, **kwargs):
                nonlocal failed
                if not failed and first_name in str(path) and path.name == "wt_reference" and path.parent.name == "repeat1":
                    failed = True
                    raise OSError("SYNTHETIC AC directory creation failure")
                return original_mkdir(path, *args, **kwargs)
            with mock.patch.object(Path, "mkdir", fail_one_directory), contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, workers=1, runner=runner)
            result = fx.load_json(output / "jobs" / first_name / "result.json")
            self.assertTrue(result["build_valid"])
            self.assertEqual(len(result["pairs"]), 5, "An exception must never leave a truncated tracked pair list")
            decision = fx.analyze(output, root)
            self.assertFalse(decision["complete"])
            self.assertEqual(decision["job_counts"], {"ok": 23, "failed": 1, "not_run": 0})
            self.assertGreater(decision["pair_counts"]["failed"], 0)
            self.assertEqual(fx.verify(output, root)["published_tables"], "PASS")
            before = len(calls)
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, workers=1, runner=runner)
            self.assertEqual(len(calls), before, "Verified failed jobs must not be silently rerun/promoted")

    def test_postflight_provenance_failure_keeps_raw_but_cannot_promote_pairs(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_MATCHED_") as tmp:
            root, output, _, binary = panel(tmp)
            fx.prepare(output, root)
            last_name = fx.read_csv(output / "input_manifest.csv")[-1]["job_id"]
            prepared = output / "jobs" / last_name / "individual_list.txt"
            original = prepared.read_bytes()
            base_runner, calls = fixture_runner()
            def runner(command, cwd, stage, timeout):
                metadata = base_runner(command, cwd, stage, timeout)
                if stage == "analyse" and last_name in str(cwd) and cwd.parent.name == "repeat5" and cwd.name == "mutant":
                    prepared.write_bytes(original + b"SYNTHETIC unexpected mutation-file change\n")
                return metadata
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, workers=1, runner=runner)
            result = fx.load_json(output / "jobs" / last_name / "result.json")
            self.assertTrue(result["build_valid"])
            self.assertEqual(len(result["pairs"]), 5)
            self.assertTrue(all(pair["status"] == "failed" for pair in result["pairs"]), "Successful numeric tables cannot override a provenance failure")
            self.assertTrue(all("gain_kcal_mol" not in pair for pair in result["pairs"]))
            with self.assertRaises(ValueError):
                fx.verify(output, root)
            prepared.write_bytes(original)
            decision = fx.analyze(output, root)
            self.assertFalse(decision["complete"])
            self.assertEqual(decision["job_counts"], {"ok": 23, "failed": 1, "not_run": 0})
            self.assertEqual(decision["pair_counts"], {"ok": 115, "failed": 5, "not_run": 0})
            self.assertEqual(fx.verify(output, root)["published_tables"], "PASS")
            self.assertTrue((output / "jobs" / last_name / "ac/repeat5/mutant/Interaction_fixture.fxout").exists())
            before = len(calls)
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(output, binary, root, workers=1, runner=runner)
            self.assertEqual(len(calls), before)


if __name__ == "__main__":
    unittest.main()
