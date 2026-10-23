"""Scientific integrity checks for the small local-regeneration experiment."""
import importlib.util
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "pipeline" / "feedback_pilot.py"


class FeedbackPilotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if SOURCE.exists():
            spec = importlib.util.spec_from_file_location("feedback_pilot", SOURCE)
            cls.pilot = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(cls.pilot)
        else:
            cls.pilot = None

    def test_missing_contacts_are_not_low_confidence(self):
        self.assertIsNotNone(self.pilot, "pilot helper implementation is missing")
        rows = [{"peptide_position": 1, "count": 0, "mean": ""},
                {"peptide_position": 2, "count": 3, "mean": .4},
                {"peptide_position": 3, "count": 1, "mean": .1},
                {"peptide_position": 4, "count": 2, "mean": .8}]
        eligible, selected = self.pilot.select_positions(rows, 2)
        self.assertEqual(eligible, [2, 3, 4])
        self.assertEqual(selected, [3, 2])

    def test_msa_removes_peptide_homologs_and_preserves_insertions(self):
        self.assertIsNotNone(self.pilot, "pilot helper implementation is missing")
        msa = "#3,2\t1,1\n>101\t102\nACDEF\n>rec\nAqCD--\n>paired\nACDEF\n>pep\n---EF\n"
        result, _ = self.pilot.make_shared_msa(msa, "ACD", "EF", "GW")
        records = self.pilot.parse_a3m(result)
        self.assertEqual(records[0][1], "ACDGW")
        self.assertIn(("rec", "AqCD--"), records)
        self.assertEqual(records[-1][1], "---GW")
        self.assertNotIn("EF", result)
        self.assertEqual(len(records), 3)

    def test_msa_rejects_wrong_parent_query(self):
        self.assertIsNotNone(self.pilot, "pilot helper implementation is missing")
        with self.assertRaises(ValueError):
            self.pilot.make_shared_msa("#3,2\t1,1\n>101\t102\nACDFF\n", "ACD", "EF", "GW")

    def test_only_selected_positions_can_change(self):
        self.assertIsNotNone(self.pilot, "pilot helper implementation is missing")
        self.pilot.validate_child("ACDEF", "AWDLF", [2, 4])
        for child in ["AWDLE", "ACDLF", "AWDLX"]:
            with self.assertRaises(ValueError):
                self.pilot.validate_child("ACDEF", child, [2, 4])


if __name__ == "__main__":
    unittest.main()
