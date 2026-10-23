"""Regression tests for scientific input integrity and the approved cloud pilot."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline"))
import step3_graphpep as gp


def atom(serial: int, name: str, residue: str, chain: str, number: int, x: float) -> str:
    return (f"ATOM  {serial:5d} {name:^4s} {residue:>3s} {chain}{number:4d}    "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}  1.00 80.00           {name[0]:>2s}\n")


class GraphPepInputsTests(unittest.TestCase):
    def test_each_complex_keeps_its_own_receptor_coordinate_frame(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            for ident, offset, resname in [("P1", 0., "ALA"), ("P2", 100., "GLY")]:
                (base / f"{ident}_complex.pdb").write_text(
                    atom(1, "CA", "ALA", "A", 1, offset) +
                    atom(2, "CB", "ALA", "A", 1, offset + .7) +
                    atom(3, "CA", resname, "B", 1, offset + 4.) + "END\n")
            manifest = pd.DataFrame([
                {"peptide_id": "P1", "peptide_sequence": "A", "source": "pilot"},
                {"peptide_id": "P2", "peptide_sequence": "G", "source": "pilot"},
            ])
            pairs = gp.build_inputs(manifest, str(base))
            self.assertIsInstance(pairs, list, "Each complex must produce its own receptor/peptide pair")
            self.assertEqual(len(pairs), 2)
            for pair, expected in zip(pairs, [0., 100.]):
                rec = Path(pair["receptor_path"]).read_text()
                pep = Path(pair["peptide_path"]).read_text()
                self.assertAlmostEqual(float(rec.splitlines()[0][30:38]), expected)
                pep_atom = next(line for line in pep.splitlines() if line.startswith("ATOM"))
                self.assertAlmostEqual(float(pep_atom[30:38]), expected + 4.)
                self.assertIn(" CB ", rec, "Sidechain atoms must survive splitting")
                self.assertEqual(pep.count("MODEL"), 1)

    def test_wrong_peptide_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "P1_complex.pdb").write_text(
                atom(1, "CA", "ALA", "A", 1, 0.) +
                atom(2, "CA", "GLY", "B", 1, 4.) + "END\n")
            manifest = pd.DataFrame([{"peptide_id": "P1", "peptide_sequence": "A", "source": "pilot"}])
            with self.assertRaisesRegex(ValueError, "sequence"):
                gp.build_inputs(manifest, str(base))


class PilotTests(unittest.TestCase):
    def test_prepare_creates_reproducible_two_chain_pilot(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            target = base / "target.fasta"
            target.write_text(">synthetic receptor\n" + "A" * 20 + "\n")
            candidates = base / "generated.csv"
            with candidates.open("w", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(["peptide_id", "peptide_sequence", "ppl_score"])
                for i, seq in enumerate(["AAAAAAAAAA", "AAAAAAAAAA", "AAAAAAAAGG", "GGGGGGGGGG", "RRRRRRRRRR", "KKKKKKKKKK"]):
                    writer.writerow([f"P{i}", seq, 1 if i == 1 else i + 1])
            out = base / "out"
            command = [sys.executable, str(ROOT / "pipeline/experiment0.py"), "prepare",
                       "--peptides", str(candidates), "--target-fasta", str(target),
                       "--receptor-range", "2-19", "--out-dir", str(out),
                       "--reference-notebook", str(ROOT / "pipeline/colabfold_batch_reference.ipynb")]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            config = json.loads((out / "pilot_config.json").read_text())
            self.assertEqual(config["n_unique_sequences"], 5)
            self.assertEqual(config["receptor_length"], 18)
            with (out / "pilot_candidates.csv").open() as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 2)
            self.assertNotEqual(rows[0]["peptide_sequence"], rows[1]["peptide_sequence"])
            fasta = (out / "pilot_complexes.fasta").read_text()
            self.assertEqual(fasta.count(":"), 2)
            notebook = json.loads((out / "Experiment0_ColabFold.ipynb").read_text())
            code = "\n".join("".join(c["source"]) for c in notebook["cells"] if c["cell_type"] == "code")
            self.assertIn("alphafold2_multimer_v3", code)
            self.assertIn(rows[0]["peptide_sequence"], code)
            self.assertIn("files.download", code)
            self.assertNotIn("drive.mount", code)
            before = (out / "pilot_candidates.csv").read_bytes()
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(before, (out / "pilot_candidates.csv").read_bytes())


class ReportTests(unittest.TestCase):
    def test_report_selects_lower_scores_and_preserves_missing_feedback(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            scores = base / "scores.csv"
            scores.write_text("peptide_id,source,binding_score,red_residue_ratio,n_contacts\n"
                              "best,ppl,-4,0,10\nworst,ppl,-1,0.5,10\n"
                              "invalid,ppl,-10,,0\nbe,bindeval,-3,,10\n")
            report = base / "report.md"
            result = subprocess.run([
                sys.executable, str(ROOT / "pipeline/step4_analyze.py"),
                "--peptides", str(base / "missing.csv"), "--scores", str(scores),
                "--top-ppl", str(base / "ppl.txt"), "--top-bindeval", str(base / "be.txt"),
                "--report", str(report), "--figure", str(base / "plot.png"), "--top-n", "1"],
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            content = report.read_text()
            self.assertIn("lower", content)
            self.assertIn("PPL filtering retained", content)
            self.assertNotIn("higher mean GraphPep", content)
            self.assertNotIn("| invalid |", content)


class ImportTests(unittest.TestCase):
    def test_swapped_chains_confidence_and_original_receptor_numbering(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            raw = base / "cloud"
            raw.mkdir()
            config = {"receptor_sequence": "AC", "receptor_range_1based_inclusive": [20, 21],
                      "prediction_method": "AFM", "candidates": [
                          {"peptide_id": "P1", "peptide_sequence": "G" * 10, "source": "pilot"}]}
            (base / "pilot_config.json").write_text(json.dumps(config))
            # The batch CLI adds a numeric job suffix; the notebook uses the bare ID.
            name = "P1_0_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_007.pdb"
            text = "".join(atom(i, "CA", "GLY", "A", i, 4. + i*.1) for i in range(1, 11))
            text += atom(11, "CA", "ALA", "B", 1, 0.) + atom(12, "CA", "CYS", "B", 2, 8.)
            (raw / name).write_text(text + "END\n")
            conf = {"plddt": [90.] * 10 + [50.] * 2,
                    "pae": [[5.] * 12 for _ in range(12)], "iptm": .7, "ptm": .8}
            (raw / name.replace("_unrelaxed_", "_scores_").replace(".pdb", ".json")).write_text(json.dumps(conf))
            result = subprocess.run([sys.executable, str(ROOT / "pipeline/experiment0.py"),
                                     "import", "--results", str(raw), "--pilot-dir", str(base)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = pd.read_csv(base / "graphpep_models.csv")
            self.assertEqual(len(manifest), 1)
            self.assertEqual(manifest.iloc[0]["peptide_plddt"], 90.)
            self.assertEqual(manifest.iloc[0]["interface_pae"], 5.)
            path = Path(manifest.iloc[0]["complex_path"])
            lines = path.read_text().splitlines()
            rec = [line for line in lines if line.startswith("ATOM") and line[21] == "A"]
            self.assertEqual([int(line[22:26]) for line in rec], [20, 21])
            pairs = gp.build_inputs(manifest, str(base / "inputs"))
            self.assertEqual(pairs[0]["receptor_map"], {1: "20:A", 2: "21:A"})
            feedback = pd.DataFrame({"DECOY_ID": [1], "RESID_PEP": ["1:B"],
                                     "RESID_PRO": ["2:A"], "CONFIDENT_SCORES": [.8]})
            mapped = gp.map_feedback(feedback, pairs[0])
            self.assertEqual(mapped.iloc[0]["RESID_PRO"], "21:A")

    def test_no_contact_stays_missing_and_does_not_run_graphpep(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            (base / "P1_complex.pdb").write_text(atom(1, "CA", "ALA", "A", 1, 0.) +
                                                atom(2, "CA", "GLY", "B", 1, 100.) + "END\n")
            manifest = pd.DataFrame([{"peptide_id": "P1", "peptide_sequence": "G", "source": "pilot"}])
            pair = gp.build_inputs(manifest, str(base))[0]
            result = gp.score_pair(pair, base / "does-not-exist", base / "missing.py", .7)
            self.assertEqual(result["graphpep_status"], "no_contact")
            self.assertTrue(pd.isna(result["binding_score"]))
            self.assertTrue(pd.isna(result["red_residue_ratio"]))


if __name__ == "__main__":
    unittest.main()
