"""Guard scientific sign, pairing and chemical-sequence mistakes."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from affinity_reference import REFERENCES, affinity_rank_correlation, paired_median_gain


def test_affinity_correlation_sign():
    assert affinity_rank_correlation([0.5, 5, 500], [-4, -3, -2]) == 1.0
    assert affinity_rank_correlation([0.5, 5, 500], [-2, -3, -4]) == -1.0


def test_paired_median_is_not_difference_of_medians():
    # A difference of medians would give +0.5 for these same observations.
    assert abs(paired_median_gain([-4, -3, -2], [-3.9, -2.5, -3.5]) + 0.1) < 1e-12


def test_reference_keeps_assayed_length_units_and_mutation_positions():
    references = {row.peptide_id: row for row in REFERENCES}
    assert len(references) == 10
    assert references["P53_WT"].sequence == "ETFSDLWKLLPE"
    assert references["P53_WT"].kd_nM == 440
    assert references["PMI_N8A"].kd_nM == 0.49
    for row in REFERENCES:
        assert len(row.sequence) == 12
        parent = references[row.parent_id]
        changes = [i for i, (left, right) in enumerate(zip(parent.sequence, row.sequence), 1) if left != right]
        assert changes == ([row.position] if row.position else [])


def test_new_structure_input_does_not_reuse_another_structures_score():
    from affinity_diagnostic import score_input_metadata
    row = {"peptide_id": "PMI_WT", "sequence": "TSFAEYWNLLSP", "seed": 20261011,
           "binding_score": -9., "contact_feedback_path": "old_full_structure.csv",
           "structure_sha256": "old_hash", "graphpep_status": "ok", "n_contacts": 100,
           "n_geometry_contact_pairs": 100, "graphpep_source_signature": "old_signature",
           "contact_pairs_6A": 100}
    cleaned = score_input_metadata(row)
    assert cleaned == {"peptide_id": "PMI_WT", "sequence": "TSFAEYWNLLSP", "seed": 20261011}


if __name__ == "__main__":
    for check in (test_affinity_correlation_sign,
                  test_paired_median_is_not_difference_of_medians,
                  test_reference_keeps_assayed_length_units_and_mutation_positions,
                  test_new_structure_input_does_not_reuse_another_structures_score):
        check()
        print(f"PASS {check.__name__}")
