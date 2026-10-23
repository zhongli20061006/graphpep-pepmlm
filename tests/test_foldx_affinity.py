"""Scientific guards; synthetic executables/energies stay in TemporaryDirectory ONLY."""
from pathlib import Path
import contextlib
import io
import math
import shutil
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))
import foldx_affinity as fx
import verify_foldx_scores as independent_score_audit


def atom(serial, name, residue, chain, number, element, altloc=" "):
    return (f"ATOM  {serial:5d} {name:^4s}{altloc}{residue:3s} {chain}{number:4d}    "
            f"{serial * .1:8.3f}{number * .2:8.3f}{-serial * .3:8.3f}{1.:6.2f}{90.:6.2f}          {element:>2s}  ")


def synthetic_pdb(expected):
    """Fixture coordinates are arbitrary and NEVER presented as molecular predictions."""
    reverse = {value: key for key, value in fx.AA.items()}
    lines, serial = [], 1
    for chain, sequence in expected.items():
        for number, letter in enumerate(sequence, 1):
            for name, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("H", "H")):
                lines.append(atom(serial, name, reverse[letter], chain, number, element))
                serial += 1
            if number == len(sequence):
                lines.append(atom(serial, "OXT", reverse[letter], chain, number, "O"))
                serial += 1
        lines.append("TER")
    return "\n".join(lines + ["END", ""])


def synthetic_panel(directory):
    source, output = directory / "SYNTHETIC_SOURCE", directory / "SYNTHETIC_OUTPUT"
    (source / "analysis").mkdir(parents=True)
    # Reference chemical identities/Kd are copied unchanged; all structures are synthetic.
    shutil.copyfile(fx.SOURCE / "reference_affinity.csv", source / "reference_affinity.csv")
    shutil.copyfile(fx.SOURCE / "protocol.json", source / "protocol.json")
    rows = fx.read_csv(fx.SOURCE / "analysis/model_results.csv")
    for row in rows:
        path = source / f"{row['peptide_id']}_seed{row['seed']}_model{row['model_number']}.pdb"
        path.write_text(synthetic_pdb({"A": fx.RECEPTOR, "B": row["sequence"]}), encoding="ascii")
        row.update(complex_path=str(path), structure_sha256=fx.sha(path))
    (source / "analysis/model_results.csv").write_text(fx.csv_text(rows), encoding="utf-8")
    return source, output


class ParserGuards(unittest.TestCase):
    def table(self, directory, text, name="Interaction_fixture.fxout"):
        path = Path(directory) / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_header_sign_prologue_and_reversed_chains(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_FOLDX_PARSE_") as directory:
            for chains, energy in (("A\tB", -9.75), ("B\tA", 4.25)):
                path = self.table(directory, f"FoldX version 5\n9999 1234 0000\nPdb\tGroup1\tGroup2\tStability\tInteraction Energy\ninput_Repair.pdb\t{chains}\t-999\t{energy}\n")
                self.assertEqual(fx.parse_interaction(path, "input_Repair.pdb"), energy)

    def test_wrong_chain_nonfinite_multiple_and_stability_rejected(self):
        header = "Pdb\tGroup1\tGroup2\tInteraction Energy\n"
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_FOLDX_PARSE_") as directory:
            invalid = [header + "input_Repair.pdb\tA\tC\t-3\n", header + "input_Repair.pdb\tAB\tB\t-3\n",
                       header + "input_Repair.pdb\tA\tB\tnan\n", header + "input_Repair.pdb\tA\tB\tinf\n",
                       header + "input_Repair.pdb\tA\tB\t-inf\n", header + "input_Repair.pdb\tA\tB\t-3\n" * 2,
                       "Pdb\tGroup1\tGroup2\tStability\ninput_Repair.pdb\tA\tB\t-3\n",
                       "Pdb\tInteraction Energy\ninput_Repair.pdb\t-3\n", header + "wrong.pdb\tA\tB\t-3\n"]
            for text in invalid:
                with self.subTest(text=text):
                    with self.assertRaises((ValueError, KeyError)):
                        fx.parse_interaction(self.table(directory, text), "input_Repair.pdb")
            with self.assertRaises(ValueError):
                fx.parse_interaction(self.table(directory, header + "input_Repair.pdb\tA\tB\t-3\n", "Dif_BuildModel.fxout"))

    def test_summary_is_accepted_only_with_named_interaction_field(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_FOLDX_PARSE_") as directory:
            path = self.table(directory, "Summary\nGroup2\tInteraction Energy\tPdb\tGroup1\nB\t-2.25\tinput_Repair.pdb\tA\n", "Summary_fixture.fxout")
            self.assertEqual(fx.parse_interaction(path, "input_Repair.pdb"), -2.25)


class ScientificGuards(unittest.TestCase):
    def test_paired_median_not_difference_of_medians(self):
        parent, variant = {1: -4., 2: -3., 3: -2.}, {1: -3.9, 2: -2.5, 3: -3.5}
        gain, matched = fx.paired_median_gain(parent, variant)
        self.assertAlmostEqual(gain, -.1)
        self.assertEqual(matched, [1, 2, 3])
        self.assertAlmostEqual(fx.statistics.median(parent.values()) - fx.statistics.median(variant.values()), .5)
        self.assertEqual(fx.paired_median_gain({1: -4.}, {1: -5., 2: -7.}), (None, [1]))

    def test_hydrogens_removed_without_changing_heavy_atoms(self):
        expected = {"A": "G", "B": "A"}
        original = synthetic_pdb(expected)
        original = original.replace("TER\nEND", atom(999, "1D", "ALA", "B", 1, "") + "\nTER\nEND")
        exported, counts = fx.heavy_pdb(original, expected)
        heavy = [line for line in original.splitlines() if line.startswith("ATOM") and line[12:16].strip() not in ("H", "1D")]
        actual = [line for line in exported.splitlines() if line.startswith("ATOM")]
        self.assertEqual(actual, heavy)
        self.assertEqual(counts, {"heavy_atoms": 10, "removed_hydrogens": 3})
        # An explicit heavy element wins over a misleading H-prefixed atom name.
        extra = atom(1000, "HG", "ALA", "B", 1, "C")
        exported, counts = fx.heavy_pdb(original + extra + "\n", expected)
        self.assertIn(extra, exported)
        self.assertEqual(counts["heavy_atoms"], 11)

    def test_duplicate_altloc_sequence_missing_backbone_terminus_rejected(self):
        expected = {"A": "G", "B": "A"}
        original = synthetic_pdb(expected)
        lines = original.splitlines()
        variants = [original + lines[0] + "\n", original.replace(lines[0], lines[0][:16] + "A" + lines[0][17:]),
                    original.replace("GLY", "UNK"), original.replace(lines[0] + "\n", ""),
                    "\n".join(line for line in lines if line[12:16].strip() != "OXT") + "\n"]
        for variant in variants:
            with self.subTest(variant=variant[:100]):
                with self.assertRaises(ValueError):
                    fx.heavy_pdb(variant, expected)


class WorkflowGuards(unittest.TestCase):
    def test_missing_binary_has_explicit_pending_state_and_no_energies(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_FOLDX_MISSING_BINARY_") as directory:
            source, output = synthetic_panel(Path(directory))
            fx.prepare(source, output)
            for name in ("repair.cfg", "analyse.cfg"):
                config = (output / "jobs/PMI_WT_seed20261011_model1" / name).read_text()
                self.assertIn("noNterm=false\n", config)
                self.assertIn("noCterm=false\n", config)
            frozen = fx.sha(output / "protocol.json")
            fx.prepare(source, output)
            self.assertEqual(frozen, fx.sha(output / "protocol.json"))
            with self.assertRaisesRegex(ValueError, "no energies produced"):
                fx.run(source, output, Path(directory) / "ABSENT_OFFICIAL_FOLDX")
            self.assertEqual(fx.load_json(output / "execution_state.json")["state"], "awaiting_foldx")
            self.assertEqual(list(output.rglob("*.fxout")), [])
            self.assertEqual(list(output.rglob("result.json")), [])
            self.assertFalse((output / "analysis").exists())
            verified = fx.verify(source, output)
            self.assertEqual(verified["real_scoring"], "NOT_RUN")
            self.assertFalse(verified["scoring_complete"])
            self.assertFalse(verified["continuation_pass"])
            # Changed source structures cannot silently reuse a frozen preparation.
            first = Path(fx.read_csv(output / "input_manifest.csv")[0]["source_path"])
            first.write_text(first.read_text() + "REMARK changed source\n")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                fx.prepare(source, output)

    @unittest.skipIf(sys.platform == "win32", "Synthetic Linux executable fixture runs in authorized WSL Python")
    def test_synthetic_run_raw_audit_failures_resume_and_tool_change(self):
        with tempfile.TemporaryDirectory(prefix="SYNTHETIC_FOLDX_WORKFLOW_") as directory:
            directory = Path(directory)
            source, output = synthetic_panel(directory)
            binary = directory / "SYNTHETIC_TEST_EXECUTABLE"
            binary.write_text(f"#!{sys.executable}\n" + '''# SYNTHETIC TEST FIXTURE, NOT FOLDX OR A MOLECULAR PREDICTION.
from pathlib import Path
import math
import sys
if '--help' in sys.argv:
    print('FoldX version 5.0 SYNTHETIC TEST FIXTURE')
    raise SystemExit(0)
assert sys.argv[1] == '-f'
config = Path(sys.argv[2]).read_text()
if 'command=RepairPDB' in config:
    if Path.cwd().name.startswith('P53_K24A_'):
        print('Intentional synthetic repair failure', file=sys.stderr)
        raise SystemExit(7)
    Path('input_Repair.pdb').write_bytes(Path('input.pdb').read_bytes())
else:
    kd = {'PMI_WT':3.2,'PMI_N8A':.49,'PMI_Y6A':610.,'PMI_F3A':38000.,'PMI_W7A':160000.,'P53_WT':440.,'P53_P27A':51.,'P53_K24A':230.,'P53_L22A':5000.,'P53_L26A':27000.}
    pid = Path.cwd().name.split('_seed')[0]
    energy = math.log10(kd[pid]) - 8. + .05 * int(Path.cwd().name.rsplit('model',1)[1])
    content = f'SYNTHETIC TEST OUTPUT ONLY\\nPdb\\tGroup1\\tGroup2\\tInteraction Energy\\tStability\\ninput_Repair.pdb\\tA\\tB\\t{energy}\\t-999.0\\n'
    Path('Interaction_fixture.fxout').write_text(content)
    Path('Summary_fixture.fxout').write_text(content)
''', encoding="ascii")
            binary.chmod(0o755)
            with contextlib.redirect_stdout(io.StringIO()):
                state = fx.run(source, output, binary)
            self.assertEqual((state["scored_models"], state["failed_models"]), (54, 6))
            report = fx.analyze(source, output)
            self.assertEqual(report["real_scoring"], "COMPLETE_WITH_FAILURES")
            self.assertTrue(report["scoring_complete"])
            self.assertFalse(report["continuation_pass"])
            self.assertAlmostEqual(report["seed_gates"][0]["rho"]["overall"], 1.)
            # Perfect correlation and seven correct/comparable main directions
            # cannot pass when the descriptive K24A row leaves only nine peptides.
            for gate in report["seed_gates"]:
                self.assertTrue(gate["all_seven_comparable"])
                self.assertEqual(gate["correct_of_seven"], 7)
                self.assertEqual(gate["n_comparable_peptides"], 9)
                self.assertFalse(gate["all_ten_eligible"])
                self.assertFalse(gate["gate_pass"])
            # Model 3 is deliberately absent: provisional seed gates may pass,
            # but a 40/60 partial run must never authorize continuation.
            partial_source, partial_output = synthetic_panel(directory / "PARTIAL_GATE_FIXTURE")
            partial_binary = directory / "SYNTHETIC_PARTIAL_EXECUTABLE"
            partial_binary.write_text(binary.read_text().replace("if Path.cwd().name.startswith('P53_K24A_'):", "if False:  # no failed jobs in partial snapshot fixture"))
            partial_binary.chmod(0o755)
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(partial_source, partial_output, partial_binary)
            partial_state = fx.load_json(partial_output / "execution_state.json")
            for record in fx.read_csv(partial_output / "input_manifest.csv"):
                if record["model_number"] == "3":
                    job = partial_output / "jobs" / record["job_id"]
                    for name in [*fx.job_artifacts(job), "result.json"]:
                        (job / name).unlink()
                    del partial_state["jobs"][record["job_id"]]
            partial_state.update(state="running", scored_models=40, failed_models=0)
            fx.save_json(partial_output / "execution_state.json", partial_state)
            partial = fx.analyze(partial_source, partial_output)
            self.assertEqual(partial["counts"], {"ok": 40, "failed": 0, "not_run": 20})
            self.assertTrue(all(gate["gate_pass"] for gate in partial["seed_gates"]))
            self.assertFalse(partial["scoring_complete"])
            self.assertFalse(partial["continuation_pass"])
            self.assertEqual(fx.verify(source, output)["published_tables"], "PASS")
            with mock.patch.object(independent_score_audit, "SOURCE", source), mock.patch.object(independent_score_audit, "OUTPUT", output), contextlib.redirect_stdout(io.StringIO()):
                independent_score_audit.main()
            audit = fx.load_json(output / "analysis/independent_score_verification.json")
            self.assertEqual((audit["status"], audit["real_scored_models"], audit["failed_models_retained"]), ("passed", 54, 6))
            self.assertFalse(audit["affinity_screen_pass"])
            job = output / "jobs/PMI_WT_seed20261011_model1"
            raw_hash = fx.sha(job / "repair.command.json")
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(source, output, binary)
            self.assertEqual(fx.sha(job / "repair.command.json"), raw_hash)
            # Separate synthetic panel exercises the quality/matched-model rules.
            quality_source, quality_output = synthetic_panel(directory / "QUALITY_SCOPE_FIXTURE")
            rows = fx.read_csv(quality_source / "analysis/model_results.csv")
            for row in rows:
                key = (row["peptide_id"], int(row["seed"]), int(row["model_number"]))
                if key in {("PMI_WT", 20261011, 1), ("PMI_N8A", 20261011, 2), ("P53_P27A", 20261012, 1), ("P53_P27A", 20261012, 2)}:
                    row.update(iptm="0.6", structure_pass="False")
            (quality_source / "analysis/model_results.csv").write_text(fx.csv_text(rows), encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                fx.run(quality_source, quality_output, binary)
            fx.analyze(quality_source, quality_output)
            directions = fx.read_csv(quality_output / "analysis/mutation_directions.csv")
            n8a = next(r for r in directions if r["scope"] == "primary_quality" and r["peptide_id"] == "PMI_N8A" and r["seed"] == "20261011")
            self.assertEqual((n8a["n_matched_models"], n8a["comparable"]), ("1", "False"))
            summaries = fx.read_csv(quality_output / "analysis/sequence_summary.csv")
            p27a = next(r for r in summaries if r["scope"] == "primary_quality" and r["peptide_id"] == "P53_P27A" and r["seed"] == "20261012")
            self.assertEqual((p27a["n_scored_selected"], p27a["eligible"]), ("1", "False"))
            all_p27a = next(r for r in summaries if r["scope"] == "all_models" and r["peptide_id"] == "P53_P27A" and r["seed"] == "20261012")
            self.assertEqual((all_p27a["n_scored_selected"], all_p27a["eligible"]), ("3", "True"))
            table = output / "analysis/mutation_directions.csv"
            table.write_text(table.read_text() + "tampered\n")
            with self.assertRaisesRegex(ValueError, "Published analysis differs"):
                fx.verify(source, output)
            fx.analyze(source, output)
            binary.write_text(binary.read_text() + "# Changed tool fingerprint\n")
            with self.assertRaisesRegex(ValueError, "fingerprint changed"):
                fx.run(source, output, binary)


if __name__ == "__main__":
    unittest.main(verbosity=2)
