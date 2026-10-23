"""Small hand-checkable cases for scientifically consequential MSA handling."""
import importlib.util
import unittest
from pathlib import Path


class MsaCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).with_name("msa_core.py")
        if not path.exists():
            cls.core = None
            return
        spec = importlib.util.spec_from_file_location("msa_core", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        cls.core = module

    def setUp(self):
        self.assertIsNotNone(self.core, "MSA implementation has not been written")

    def test_projection_excludes_tail_and_chain_boundary_insertions(self):
        # Full reference A C D E F; subject A W C D Y E F.
        # Core is ACD: W is inside core, Y is at its right boundary.
        projected = self.core.project_core("ACDEF", "AWCDYEF", [[0, 1, 1, 3, 3, 5], [0, 1, 2, 4, 5, 7]], 3)
        self.assertEqual(projected["a3m"], "AwCD")
        self.assertEqual(projected["match"], "ACD")
        self.assertEqual(projected["source_positions"], [1, 3, 4])

    def test_projection_records_deletion_without_shifting_numbering(self):
        projected = self.core.project_core("ACDEF", "ACEF", [[0, 2, 3, 5], [0, 2, 2, 4]], 4)
        self.assertEqual(projected["a3m"], "AC-E")
        self.assertEqual(projected["source_positions"], [1, 2, None, 3])

    def test_identity_does_not_count_gap_gap(self):
        self.assertEqual(self.core.identity("AC--", "AC--"), 0.5)
        self.assertEqual(self.core.identity("ACDE", "ACDF"), 0.75)

    def test_neff_counts_redundancy_and_greedy_clusters(self):
        result = self.core.diversity(["AAAAA", "AAAAC", "CCCCC"])
        self.assertEqual(result["greedy_clusters_80"], 2)
        self.assertAlmostEqual(result["neff_80"], 2.0)
        self.assertEqual(result["neighbor_counts_80"], [2, 2, 1])
        self.assertEqual(result["column_nongap_fraction"], [1.0] * 5)

    def test_duplicates_cannot_inflate_effective_depth(self):
        self.assertAlmostEqual(self.core.diversity(["AAAAA"] * 20)["neff_80"], 1.0)

    def test_farthest_first_selects_diversity_before_near_duplicate(self):
        records = [dict(accession="B", match="AAAAC"), dict(accession="A", match="CCCCC")]
        chosen = self.core.farthest_first("AAAAA", records, 1)
        self.assertEqual([r["accession"] for r in chosen], ["A"])

    def test_complex_unpaired_rows_have_correct_chain_lengths(self):
        content = self.core.complex_a3m("ACD", "EF", [("homolog", "AwC-")])
        lines = content.splitlines()
        self.assertEqual(lines[0], "#3,2\t1,1")
        rows = [lines[i + 1] for i in range(1, len(lines), 2)]
        self.assertEqual(rows, ["ACDEF", "ACD--", "AwC---", "---EF"])
        self.assertTrue(all(sum(not c.islower() for c in row) == 5 for row in rows))

    def test_readthrough_and_fragments_are_excluded(self):
        base = {"Fragment": "", "Protein names": "Bcl-2-like protein 2", "Gene Names": "BCL2L2", "Sequence": "A" * 193, "Length": "193"}
        self.assertIsNone(self.core.metadata_reason(base))
        self.assertEqual(self.core.metadata_reason(dict(base, Fragment="fragment")), "EXCLUDE_FRAGMENT")
        self.assertEqual(self.core.metadata_reason(dict(base, **{"Gene Names": "BCL2L2-PABPN1"})), "EXCLUDE_FUSION_READTHROUGH")

    def test_floor_failure_is_retained(self):
        floor = {"unique_nonquery_homolog_rows_min": 8, "greedy_clusters_at_80pct_min": 3, "neff_at_80pct_min": 3.0, "each_core_column_nongap_fraction_min": 0.8}
        passed = self.core.readiness(20, self.core.diversity(["AAAAA"] * 21), floor)
        self.assertFalse(all(passed.values()))
        self.assertFalse(passed["greedy_clusters"])
        self.assertFalse(passed["neff"])


if __name__ == "__main__":
    unittest.main()
