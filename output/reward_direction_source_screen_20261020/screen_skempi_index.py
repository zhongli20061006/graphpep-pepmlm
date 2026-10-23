"""Create a coarse SKEMPI v2.0 index screen; no labels are released here."""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "raw" / "skempi_v2.csv"
OUTPUT = ROOT / "skempi_index_screen.json"


def point_value(raw: str | None) -> float | None:
    """Return finite exact positive affinity values; retain censored rows as missing."""
    if raw is None or not raw.strip():
        return None
    value = raw.strip()
    if any(mark in value for mark in (">", "<", "~", "≈")):
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) and number > 0 else None


def main() -> None:
    with SOURCE.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter=";"))

    method_counts = Counter(row.get("Method", "") for row in rows)
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (
            row.get("#Pdb", ""),
            row.get("Reference", ""),
            row.get("Method", ""),
            row.get("Protein 1", ""),
            row.get("Protein 2", ""),
        )
        grouped[key].append(row)

    systems: list[dict[str, object]] = []
    peptide_terms = ("peptide", "peptid", "sh3", "pdz", "evh1", "ww domain", "motif")
    keyword_systems: list[dict[str, object]] = []
    for key, records in grouped.items():
        gains = 0
        losses = 0
        finite = 0
        variants: list[str] = []
        mutation_fields: list[str] = []
        for row in records:
            mut = point_value(row.get("Affinity_mut (M)"))
            wt = point_value(row.get("Affinity_wt (M)"))
            if mut is None or wt is None:
                continue
            ratio = mut / wt
            finite += 1
            if ratio <= 0.5:
                gains += 1
            elif ratio >= 2.0:
                losses += 1
            variants.append(row.get("Mutation(s)_cleaned", ""))
            mutation_fields.append(row.get("Mutation(s)_PDB", ""))
        record = {
            "pdb_chain_pair": key[0],
            "reference_id": key[1],
            "method_field": key[2],
            "protein_1": key[3],
            "protein_2": key[4],
            "rows": len(records),
            "finite_exact_affinity_rows": finite,
            "point_ratio_gain_le_0_5": gains,
            "point_ratio_loss_ge_2": losses,
            "both_direction_point_candidates": gains > 0 and losses > 0,
            "mutations_cleaned": sorted(set(variants)),
            "mutations_pdb_notation": sorted(set(mutation_fields)),
            "label_status": "index_only_unverified_against_primary_paper",
        }
        systems.append(record)
        names = f"{key[3]} {key[4]}".lower()
        if any(term in names for term in peptide_terms):
            keyword_systems.append(record)

    systems.sort(
        key=lambda item: (
            bool(item["both_direction_point_candidates"]),
            int(item["finite_exact_affinity_rows"]),
            int(item["rows"]),
        ),
        reverse=True,
    )
    keyword_systems.sort(
        key=lambda item: (
            bool(item["both_direction_point_candidates"]),
            int(item["finite_exact_affinity_rows"]),
            int(item["rows"]),
        ),
        reverse=True,
    )
    result = {
        "source": "SKEMPI v2.0 official CSV",
        "source_url": "https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv",
        "screen_scope": "coarse index only; source table, assay, chain, peptide chemistry, uncertainty and censoring semantics not yet adjudicated",
        "raw_rows": len(rows),
        "group_key": ["PDB chain pair", "reference id", "method field", "protein 1", "protein 2"],
        "point_ratio": "Affinity_mut (M) / Affinity_wt (M); only raw finite exact values; no censored/imputed rows",
        "point_cutoffs_used_for_lead_screen_only": {"gain": "<=0.5", "loss": ">=2.0"},
        "unique_group_count": len(systems),
        "groups_with_both_point_directions": sum(bool(item["both_direction_point_candidates"]) for item in systems),
        "both_direction_groups": [item for item in systems if item["both_direction_point_candidates"]],
        "keyword_groups_count": len(keyword_systems),
        "method_field_counts": dict(method_counts.most_common()),
        "top_groups": systems[:100],
        "keyword_groups": keyword_systems,
        "warning": "This ranking is not a candidate qualification. SKEMPI mixes protein-protein mutation entries and is not an assay-homogeneous peptide Kd benchmark. Every lead needs chain-level PDB and primary-reference review.",
    }
    OUTPUT.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "raw_rows", "unique_group_count", "groups_with_both_point_directions",
        "keyword_groups_count", "method_field_counts"
    )}, indent=2, ensure_ascii=False))
    for item in keyword_systems[:80]:
        print(json.dumps(item, ensure_ascii=False))


if __name__ == "__main__":
    main()
