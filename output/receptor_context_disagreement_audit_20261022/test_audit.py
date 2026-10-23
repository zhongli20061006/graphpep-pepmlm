from __future__ import annotations
import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
if importlib.util.find_spec("audit_core"):
    import audit_core as core
else:
    core = None


class AuditTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(core, "audit_core must implement the frozen CPU audit")

    def test_directional_opposition_is_distinct_from_abstention(self):
        self.assertEqual(core.disagreement_type("gain", "loss"), "strict_opposite")
        self.assertEqual(core.disagreement_type("loss", "abstain"), "H_abstains")
        self.assertEqual(core.disagreement_type("tie", "gain"), "d_tie_H_directional")
        self.assertEqual(core.disagreement_type("gain", "gain"), "same_direction")

    def test_exact_pair_requires_both_parent_and_variant(self):
        model = {"parent_sequence": "TSFAEYWNLLSP", "variant_sequence": "TSFAEYWALLSP"}
        experiment = {"parent_sequence": "TSFAEYWNLASP", "variant_sequence": "TSFAEYWALASP"}
        self.assertFalse(core.exact_pair(model, experiment))
        self.assertTrue(core.exact_pair(model, dict(model)))

    def test_p53_13mer_cannot_inherit_12mer_label(self):
        model = {"parent_sequence": "ETFSDLWKLLPEN", "variant_sequence": "ETFSDAWKLLPEN"}
        experiment = {"parent_sequence": "ETFSDLWKLLPE", "variant_sequence": "ETFSDAWKLLPE"}
        self.assertFalse(core.exact_pair(model, experiment))
        self.assertEqual(core.near_pair_reason(model, experiment), "HOLD_P53_13MER_VS_12MER")

    def test_shared_receptor_name_cannot_override_full_sequence(self):
        a = core.distribution_key("weight", "RR", "PEPTIDE", 3)
        b = core.distribution_key("weight", "RS", "PEPTIDE", 3)
        self.assertNotEqual(a, b)

    def test_single_mutation_sequence_is_an_identity_key(self):
        self.assertEqual(core.variant_sequence("ACDE", 2, "G"), "AGDE")
        with self.assertRaises(ValueError):
            core.variant_sequence("ACDE", 0, "G")

    def test_no_near_match_between_l10a_and_wt_parent(self):
        a = {"parent_sequence": "TSFAEYWNLLSP", "variant_sequence": "TSFAEYWALLSP"}
        b = {"parent_sequence": "TSFAEYWNLASP", "variant_sequence": "TSFAEYWALASP"}
        self.assertEqual(core.near_pair_reason(a, b), "HOLD_PMI_L10A_BACKGROUND_MISMATCH")

    def test_l10a_lookup_covers_complete_existing_table_not_only_curated_examples(self):
        import audit_disagreements as audit
        rows = audit.read_csv(audit.SOURCES / "literature_table2.csv")
        expected = {row["sequence"] for row in rows if row["sequence"] != audit.PMI_L10A}
        experiments = audit.experimental_inventory("MDM2_FULL_SEQUENCE")
        actual = {row["variant_sequence"] for row in experiments if row["source"] == "Li2021"}
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main(verbosity=2)
