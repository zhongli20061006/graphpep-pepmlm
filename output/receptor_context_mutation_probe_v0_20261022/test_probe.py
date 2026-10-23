"""Model-free regressions for the authorized context-contrast probe."""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys
import tempfile
import time
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
if importlib.util.find_spec("probe_core"):
    import probe_core as core
else:
    core = None


class ProbeTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(core, "probe_core must implement the frozen experiment")

    def test_single_mask_has_correct_boundary_and_preserves_receptor(self):
        original = [0, 4, 5, 6, 7, 8, 2]
        masked, index = core.mask_tokens(original, [4, 5, 6], [7, 8], 1, 0, 2, 32)
        self.assertEqual(index, 4)
        self.assertEqual(masked, [0, 4, 5, 6, 32, 8, 2])
        self.assertEqual(original, [0, 4, 5, 6, 7, 8, 2])
        last, index = core.mask_tokens(original, [4, 5, 6], [7, 8], 2, 0, 2, 32)
        self.assertEqual(index, 5)
        self.assertEqual(last[-2:], [32, 2])

    def test_bad_sequence_and_site_are_rejected(self):
        for site in (0, 3):
            with self.assertRaises(ValueError):
                core.mask_tokens([0, 4, 7, 8, 2], [4], [7, 8], site, 0, 2, 32)
        with self.assertRaises(ValueError):
            core.mask_tokens([0, 4, 7, 9, 2], [4], [7, 8], 1, 0, 2, 32)

    def test_log_odds_formula_and_self_substitution(self):
        values = core.contrast({"A": -2., "G": -1.}, {"A": -2., "G": -2.5},
                               {"A": -2., "G": -1.5}, "A", "G")
        self.assertAlmostEqual(values["d_target"], 1.)
        self.assertAlmostEqual(values["h_control_1"], 1.5)
        self.assertAlmostEqual(values["h_control_2"], .5)
        self.assertAlmostEqual(values["h"], 1.)
        self.assertTrue(values["passes_filter"])
        zero = core.contrast({"A": -2.}, {"A": -3.}, {"A": -4.}, "A", "A")
        self.assertEqual(zero["h"], 0.)
        self.assertFalse(zero["passes_filter"])

    def test_control_disagreement_abstains_and_nan_fails(self):
        row = core.contrast({"A": 0., "G": 1.}, {"A": 0., "G": 2.},
                            {"A": 0., "G": -1.}, "A", "G")
        self.assertTrue(row["control_sign_disagreement"])
        self.assertFalse(row["passes_filter"])
        with self.assertRaises(ValueError):
            core.contrast({"A": 0., "G": math.nan}, {"A": 0., "G": 1.},
                          {"A": 0., "G": 1.}, "A", "G")

    def test_suggestion_limits_and_historical_exclusions(self):
        rows = [dict(parent_id="P", site=i, alt=aa, h=float(20-i), passes_filter=True,
                     excluded_historical=False) for i in range(1, 10) for aa in "AG"]
        rows[0]["excluded_historical"] = True
        chosen = core.select_suggestions(rows)
        self.assertEqual(len(chosen), 6)
        self.assertEqual(len({r["site"] for r in chosen}), 6)
        self.assertEqual(chosen[0]["alt"], "G")

    def test_job_budget_is_45_plus_15_unique_job_ids(self):
        inputs = [{"parent_id": "P", "eligible_positions": list(range(1, 9))},
                  {"parent_id": "Q", "eligible_positions": list(range(1, 8))}]
        jobs = core.build_jobs(inputs, ["target", "shuffle1", "shuffle2"])
        self.assertEqual(len(jobs), 60)
        self.assertEqual(len({j["job_id"] for j in jobs}), 60)
        self.assertEqual(sum(j["phase"] == "base" for j in jobs), 45)
        self.assertEqual([(j["parent_id"], j["site"]) for j in jobs[-15:]],
                         [(j["parent_id"], j["site"]) for j in jobs[:45]
                          if j["context_id"] == "target"][::-1])
        with self.assertRaises(ValueError):
            core.build_jobs(inputs, ["target"] * 3)

    def test_watchdog_times_once_and_stops_only_its_child(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            started = time.monotonic()
            result = core.run_bounded([sys.executable, "-c", "import time; time.sleep(.2)"],
                                      root, 2., 100000, root / "normal.log")
            actual = time.monotonic() - started
            self.assertEqual(result["status"], "completed")
            self.assertLessEqual(result["wall_seconds"], actual + .01)
            self.assertGreater(result["wall_seconds"], .18)
            result = core.run_bounded([sys.executable, "-c", "import time; time.sleep(3)"],
                                      root, .2, 100000, root / "timeout.log")
            self.assertEqual(result["status"], "wall_cap")
            self.assertLess(result["wall_seconds"], 1.2)

    def test_watchdog_detects_output_cap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            code = "from pathlib import Path; import time; Path('large').write_bytes(b'x'*4096); time.sleep(3)"
            result = core.run_bounded([sys.executable, "-c", code], root, 2., 1024,
                                      root / "disk.log")
            self.assertEqual(result["status"], "disk_cap")


if __name__ == "__main__":
    unittest.main(verbosity=2)
