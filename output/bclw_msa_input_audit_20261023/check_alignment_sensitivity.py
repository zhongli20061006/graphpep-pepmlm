"""Bounded CPU sensitivity audit; original selection/A3M remain immutable."""
import csv
import json
import time
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

import msa_core

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]


def main() -> None:
    started = time.monotonic()
    scope = json.loads((OUT / "alignment_sensitivity_scope.json").read_text())
    reference = json.loads((ROOT / "output/structure_input_audit_bclw_20261020/sources/Q92843_current.json").read_text())["sequence"]["value"]
    selected = list(csv.DictReader((OUT / "selected_homologs.csv").open()))
    rows = [reference[:164]] + [row["match"] for row in selected]
    floor = json.loads((OUT / "diversity_and_readiness.json").read_text())["frozen_floor"]
    aligner = Align.PairwiseAligner(mode="global")
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.epsilon = 1e-6
    evidence, floor_failures = [], []
    evaluated = 0
    for i, row in enumerate(selected, 1):
        alignments = aligner.align(reference, row["full_sequence"])
        total = len(alignments)
        comparisons = []
        core_rows, a3m_rows, positions = set(), set(), set()
        for j in range(min(total, scope["optimal_alignments_per_record_cap"])):
            projected = msa_core.project_core(reference, row["full_sequence"], alignments[j].coordinates.tolist(), 164)
            core_rows.add(projected["match"])
            a3m_rows.add(projected["a3m"])
            positions.update(k + 1 for k, (a, b) in enumerate(zip(projected["match"], row["match"])) if a != b)
            candidate_rows = list(rows)
            candidate_rows[i] = projected["match"]
            values = msa_core.diversity(candidate_rows)
            gates = msa_core.readiness(len(selected), values, floor)
            if not all(gates.values()):
                floor_failures.append({"accession": row["accession"], "optimum_index": j, "gates": gates})
            comparisons.append({"optimum_index": j, "match": projected["match"], "a3m": projected["a3m"],
                "greedy_clusters_80": values["greedy_clusters_80"], "neff_80": values["neff_80"],
                "minimum_column_nongap_fraction": values["minimum_column_nongap_fraction"], "gates": gates})
            evaluated += 1
        evidence.append({"accession": row["accession"], "optimal_alignment_count": total,
            "exhausted_available_optima": total <= scope["optimal_alignments_per_record_cap"],
            "distinct_projected_match_rows": len(core_rows), "distinct_projected_a3m_rows": len(a3m_rows),
            "affected_core_positions": sorted(positions), "comparisons": comparisons})
    result = {"scope": scope, "selected_homologs": len(selected),
        "selected_with_multiple_optima": sum(row["optimal_alignment_count"] > 1 for row in evidence),
        "selected_with_alternate_core_rows": sum(row["distinct_projected_match_rows"] > 1 for row in evidence),
        "selected_with_alternate_a3m_rows": sum(row["distinct_projected_a3m_rows"] > 1 for row in evidence),
        "alignment_projections_evaluated": evaluated,
        "all_available_selected_optima_enumerated": all(row["exhausted_available_optima"] for row in evidence),
        "one_row_replacement_floor_failures": floor_failures,
        "additional_execution_hold": bool(floor_failures),
        "worker_wall_seconds": time.monotonic() - started, "evidence": evidence}
    with (OUT / "alignment_sensitivity.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    print(json.dumps({key: value for key, value in result.items() if key not in {"scope", "evidence"}}, ensure_ascii=False))


if __name__ == "__main__":
    main()
