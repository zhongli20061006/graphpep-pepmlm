"""Recover source rows and source-defined BH3 comparisons without model scoring."""

from __future__ import annotations

import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent

# Visually transcribed from Table 1, main PDF p. 3 / printed p. 1963.
# The source graphic and a rendered enlargement are retained for review.
SEQUENCES = {
    "MS1": "RPEIWMTQGLRRLGDEINAYYAR",
    "MS2": "RPEIWLTQSLQRLGDEINAYYAR",
    "MS3": "RPEIWLTQHLQRLGDEINAYYAR",
    "A12": "RPEIWMGQGLRRLGDEINAYYAR",
    "B3": "RPEIWLGQSLQRLGDEINAYYAR",
    "G9": "RPEIWLGQHLQRLGDEINAYYAR",
    "NoxaA": "AELPPEFAAQLRKIGDKVYC",
    "Bim": "RPEIWIAQELRRIGDEFNAYYAR",
    "Bim_A2eT": "RPEIWITQELRRIGDEFNAYYAR",
    "Bim_A2eT_I2dM": "RPEIWMTQELRRIGDEFNAYYAR",
    "Bim_A2eT_E2gG": "RPEIWITQGLRRIGDEFNAYYAR",
    "Bim_A2eT_I3dL": "RPEIWITQELRRLGDEFNAYYAR",
    "Bim_A2eT_F4aI": "RPEIWITQELRRIGDEINAYYAR",
}

# Direction of comparison is explicitly discussed in the paper. Do not add
# reverse comparisons or select a new parent after viewing scorer outputs.
SOURCE_EDGES = [
    ("Bim", "Bim_A2eT"),
    ("Bim_A2eT", "Bim_A2eT_I2dM"),
    ("Bim_A2eT", "Bim_A2eT_E2gG"),
    ("Bim_A2eT", "Bim_A2eT_I3dL"),
    ("Bim_A2eT", "Bim_A2eT_F4aI"),
    ("A12", "MS1"),
    ("B3", "MS2"),
]


def write_json(name: str, value: object) -> None:
    (ROOT / name).write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())


def canonical_label(label: str) -> str:
    compact = re.sub(r"[\s_]+", "", label)
    if compact == "Bim23-mer":
        return "Bim"
    matches = [key for key in SEQUENCES if key.replace("_", "") == compact]
    if len(matches) != 1:
        raise ValueError(f"Unrecognized source label: {label!r}")
    return matches[0]


def parse_measurement(raw: str) -> dict[str, object]:
    if raw == "ND":
        return {"raw": raw, "status": "not_determined"}
    if raw.startswith(("<", ">")):
        return {"raw": raw, "status": "censored", "relation": raw[0], "bound_nM": float(raw[1:])}
    parts = raw.split("±")
    if len(parts) != 2:
        raise ValueError(f"Unexpected measurement: {raw!r}")
    return {
        "raw": raw,
        "status": "finite",
        "kd_nM": float(parts[0]),
        "reported_95CI_halfwidth_nM": float(parts[1]),
        "uncertainty_type": "Matlab confint on nonlinear fit to mean anisotropy curve; not replicate SD or SEM",
    }


def main() -> None:
    article = ET.parse(ROOT / "raw" / "PMC4168798.xml").getroot()
    table = article.find('.//table-wrap[@id="tbl2"]')
    if table is None:
        raise ValueError("Source Table 2 is missing")
    receptors = ["Mcl-1", "Bfl-1", "Bcl-w", "Bcl-xL", "Bcl-2", "KSBcl-2"]
    records: list[dict[str, object]] = []
    lookup: dict[tuple[str, str], dict[str, object]] = {}
    for row in table.iter("tr"):
        cells = [text(cell) for cell in row]
        if len(cells) != 7 or cells[0] == "peptide":
            continue
        peptide = canonical_label(cells[0])
        for receptor, raw in zip(receptors, cells[1:], strict=True):
            record = {"peptide": peptide, "receptor": receptor, **parse_measurement(raw)}
            records.append(record)
            lookup[(receptor, peptide)] = record
    assert len(records) == 72

    sequence_records = []
    for name, sequence in SEQUENCES.items():
        assert set(sequence) <= set("ACDEFGHIKLMNPQRSTVWY")
        assert len(sequence) == (20 if name == "NoxaA" else 23)
        sequence_records.append({
            "name": name,
            "sequence": sequence,
            "length": len(sequence),
            "source": "Table 1, visually checked in figures/Table1_sequences.png",
            "chemistry": "Bim variants: N-terminal 5/6-fluorescein amidite, C-terminal amide; NoxaA uses a different FITC-AHA construct and is not used in mutation pairs",
        })

    edge_records = []
    for parent, variant in SOURCE_EDGES:
        differences = [
            {"local_position_1based": i, "from": a, "to": b}
            for i, (a, b) in enumerate(zip(SEQUENCES[parent], SEQUENCES[variant], strict=True), 1)
            if a != b
        ]
        assert len(differences) == 1, (parent, variant, differences)
        edge_records.append({"parent": parent, "variant": variant, "mutation": differences[0]})

    pair_records = []
    counts = {}
    for receptor in ("Bcl-w", "Bcl-xL"):
        receptor_counts: Counter[str] = Counter()
        for edge in edge_records:
            parent = lookup[(receptor, str(edge["parent"]))]
            variant = lookup[(receptor, str(edge["variant"]))]
            assert parent["status"] == variant["status"] == "finite"
            w, m = float(parent["kd_nM"]), float(variant["kd_nM"])
            hw = float(parent["reported_95CI_halfwidth_nM"])
            hm = float(variant["reported_95CI_halfwidth_nM"])
            ratio = m / w
            direction = "gain" if ratio <= 0.5 else "loss" if ratio >= 2 else "below_effect_cutoff"
            independent_radius = math.hypot(hm / m, hw / w)
            adverse_radius = hm / m + hw / w
            independent = [ratio * math.exp(-independent_radius), ratio * math.exp(independent_radius)]
            adverse = [ratio * math.exp(-adverse_radius), ratio * math.exp(adverse_radius)]
            # Marginal CI endpoints are also retained as a descriptive envelope.
            # With only two marginal 95% CIs this box is not a joint 95% CI.
            endpoint_envelope = [(m - hm) / (w + hw), (m + hm) / (w - hw)]
            supports_direction = (
                direction == "gain" and max(independent[1], adverse[1]) < 1
            ) or (
                direction == "loss" and min(independent[0], adverse[0]) > 1
            )
            record = {
                "receptor": receptor,
                **edge,
                "parent_measurement": parent,
                "variant_measurement": variant,
                "ratio_variant_over_parent": ratio,
                "point_effect_direction": direction,
                "fit_conditional_approximate_log_ratio_interval_independent": independent,
                "fit_conditional_approximate_log_ratio_interval_adverse_covariance": adverse,
                "marginal_95CI_endpoint_envelope_NOT_joint95CI": endpoint_envelope,
                "direction_supported_under_declared_fit_sensitivity": supports_direction,
                "formal_benchmark_label_released": False,
            }
            pair_records.append(record)
            receptor_counts[direction] += 1
        counts[receptor] = dict(receptor_counts)

    write_json("raw_table2_measurements.json", {
        "source": "Foight et al. 2014 DOI 10.1021/cb500340w, Table 2",
        "unit": "nM",
        "cell_count": len(records),
        "status_counts": dict(Counter(str(row["status"]) for row in records)),
        "records": records,
    })
    write_json("peptide_sequences_source.json", {"records": sequence_records, "edges": edge_records})
    write_json("paired_direction_candidates.json", {
        "scope": "seven source-defined single-substitution edges, audited separately for Bcl-w and Bcl-xL; no pooled metric, no reverse edges",
        "counts": counts,
        "pairs": pair_records,
        "scores_seen": False,
    })
    write_json("fit_uncertainty_sensitivity_v1.json", {
        "version": "bh3_fit_summary_sensitivity_v1",
        "frozen_before_new_prediction_or_score": True,
        "source_semantics": "Kd fit to averaged anisotropy data from at least three replicates across three days; pointwise SD; Table 2 uncertainty is Matlab confint on the fit",
        "approximation": "Treat reported symmetric 95% CI halfwidth h as 1.96 times local fit SE. Delta-method log-ratio radius is sqrt((h_mut/K_mut)^2+(h_parent/K_parent)^2), or their sum under adverse covariance. This is conditional model-based sensitivity, not raw-replicate or biological inference.",
        "no_division_by_sqrt_n": True,
        "limits": [
            "Replicate-level Kd or raw anisotropy values and full covariance are unavailable",
            "Source fits condition on a fixed lower baseline and a peptide concentration estimated from other plateaued curves",
            "Receptor active fraction and source Bcl-w oligomerization are not included",
            "Marginal-endpoint envelopes are not claimed to have joint 95% coverage",
            "No formal benchmark labels or pilot execution are released by this sensitivity calculation",
        ],
    })

    fasta = "".join(f">{name}|CORE_ONLY_NOT_JOB_INPUT\n{sequence}\n" for name, sequence in SEQUENCES.items())
    (ROOT / "peptide_core_sequences_NOT_JOB_INPUT.fasta").write_text(fasta, encoding="utf-8")
    core_records = {}
    for name, filename, length in (
        ("Bcl-w", "Q92843_BCLW.fasta", 164),
        ("Bcl-xL", "Q07817_BCLXL.fasta", 209),
    ):
        full = "".join(line.strip() for line in (ROOT / "raw" / filename).read_text(encoding="utf-8").splitlines() if not line.startswith(">"))
        core_records[name] = {
            "native_core": full[:length], "core_length": length,
            "assay_construct": "N-terminal c-myc-tagged core after MBP/TEV cleavage; complete junction/scar not reconstructed",
            "ready_as_exact_assay_input": False,
        }
    write_json("receptor_core_audit.json", core_records)
    print(json.dumps({"cells": len(records), "sequences": len(SEQUENCES), "source_edges": len(SOURCE_EDGES), "separate_receptor_counts": counts}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
