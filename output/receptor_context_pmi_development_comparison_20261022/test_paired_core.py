"""Model-free scientific regressions for the finite PMI development comparison."""
from __future__ import annotations
import importlib.util
import math
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
if importlib.util.find_spec("paired_core"):
    import paired_core as core
else:
    core = None


class PairedTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(core, "paired_core implements the approved comparison")

    def test_kd_shared_exponent_converts_mean_and_sd(self):
        self.assertEqual(core.parse_kd("3.2±1.1E-09"), (3.2, 1.1))
        self.assertEqual(core.parse_kd("4.9±2.1E-10"), (.49, .21))
        self.assertIsNone(core.parse_kd("n.d.a"))

    def test_table_colspan_keeps_kd_columns_and_truncations(self):
        html = '<table><tr><td colspan="2">A</td><td>B<sup>1</sup></td></tr></table>'
        self.assertEqual(core.parse_tables(html), [[["A", "", "B1"]]])

    def test_rounded_printed_ratio_does_not_make_gain(self):
        self.assertEqual(core.point_direction(230., 440.), "ambiguous")
        self.assertEqual(core.point_direction(240., 440.), "ambiguous")
        self.assertEqual(core.point_direction(210., 440.), "gain")
        self.assertEqual(core.point_direction(880., 440.), "loss")

    def test_sensitivity_exposes_n8a_correlation_failure(self):
        result = core.sensitivity(.49, .21, 3.2, 1.1)
        self.assertLess(result["independent_upper"], 1.)
        self.assertGreater(result["worst_correlation_upper"], 1.)
        self.assertAlmostEqual(result["ratio"], .153125)

    def test_h_control_flip_and_tie_abstain(self):
        self.assertEqual(core.predict_d(0.), "tie")
        self.assertEqual(core.predict_d(-.01), "loss")
        self.assertEqual(core.predict_h(1., -1.), "abstain")
        self.assertEqual(core.predict_h(1., 0.), "abstain")
        self.assertEqual(core.predict_h(-1., -.2), "loss")

    def test_failures_remain_in_paired_denominator(self):
        rows = [{"point_direction": "gain", "d_prediction": "gain", "h_prediction": "abstain"},
                {"point_direction": "loss", "d_prediction": "loss", "h_prediction": "loss"},
                {"point_direction": "loss", "d_prediction": "failure", "h_prediction": "failure"}]
        result = core.paired_summary(rows)
        self.assertEqual(result["n"], 3)
        self.assertEqual(result["d_correct"], 2)
        self.assertEqual(result["h_correct"], 1)
        self.assertEqual(result["d_only_correct"], 1)
        self.assertEqual(result["neither_correct"], 1)

    def test_full_cache_key_includes_receptor_peptide_and_site(self):
        key = core.cache_key("checkpoint", "RR", "PEPTIDE", 2, "FP32_eager_TF32off")
        self.assertNotEqual(key, core.cache_key("checkpoint", "Rr", "PEPTIDE", 2, "FP32_eager_TF32off"))
        self.assertNotEqual(key, core.cache_key("checkpoint", "RR", "PEPTIDEN", 2, "FP32_eager_TF32off"))
        self.assertNotEqual(key, core.cache_key("checkpoint", "RR", "PEPTIDE", 3, "FP32_eager_TF32off"))

    def test_cached_original_table_audit_handles_blank_separator_rows(self):
        import run_development
        try:
            audit, pairs = run_development.audit_source()
        except ValueError as exc:
            self.fail(f"Original source audit must handle separator rows: {exc}")
        self.assertEqual(len(pairs), 7)
        self.assertEqual(sum(row["decision"] == "NONDATA_SEPARATOR" for row in audit), 2)
        self.assertEqual(sum(row["decision"] == "HOLD_REPLICATE_COUNT_NOT_ATTRIBUTABLE_TO_P53"
                             for row in audit), 4)
        self.assertEqual(sum(row["point_direction"] == "gain" for row in pairs), 1)
        self.assertTrue(all(row["family"] == "PMI" for row in pairs))

    def test_frozen_code_inventory_covers_imported_execution_helpers(self):
        import inspect
        import run_development
        function = getattr(run_development, "code_paths", None)
        self.assertTrue(callable(function), "Freeze must inventory imported execution dependencies")
        frozen = {path.resolve() for path in function()}
        self.assertIn(Path(inspect.getsourcefile(run_development.run_bounded)).resolve(), frozen)
        self.assertIn(Path(inspect.getsourcefile(run_development.check_manifest)).resolve(), frozen)


if __name__ == "__main__":
    unittest.main(verbosity=2)
