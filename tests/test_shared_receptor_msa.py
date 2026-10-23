"""MSA parity across parents with different peptide lengths."""
import importlib.util
from pathlib import Path
import sys
import unittest

PIPELINE = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE))


class SharedReceptorMSATests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = PIPELINE / "shared_receptor_msa.py"
        cls.helper = None
        if source.exists():
            spec = importlib.util.spec_from_file_location("shared_receptor_msa", source)
            cls.helper = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(cls.helper)

    def test_different_peptide_lengths_preserve_receptor_alignment(self):
        self.assertIsNotNone(self.helper, "variable-length shared MSA helper is missing")
        msa = "#3,2\t1,1\n>101\t102\nACDEF\n>rec\nAqCD--\n>paired\nACDEF\n>pep\n---EF\n"
        hashes = []
        for child in ["G", "GW", "GWY"]:
            with self.subTest(child=child):
                result, stats = self.helper.make_shared_msa(msa, "ACD", "EF", child)
                lines = result.splitlines()
                self.assertEqual(lines[0], f"#3,{len(child)}\t1,1")
                self.assertEqual(lines[2], "ACD" + child)
                self.assertIn("AqCD" + "-" * len(child), lines)
                self.assertEqual(lines[-1], "---" + child)
                self.assertNotIn(">paired", result)
                self.assertNotIn("EF", result)
                hashes.append(stats["receptor_alignment_sha256"])
        self.assertEqual(len(set(hashes)), 1)

    def test_rejects_changed_source_query(self):
        self.assertIsNotNone(self.helper, "variable-length shared MSA helper is missing")
        with self.assertRaises(ValueError):
            self.helper.make_shared_msa("#3,2\t1,1\n>101\t102\nACDFF\n", "ACD", "EF", "GWY")


if __name__ == "__main__":
    unittest.main()
