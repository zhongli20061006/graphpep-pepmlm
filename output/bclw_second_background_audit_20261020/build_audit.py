"""Build a source-only mutation/construct audit from archived primary records."""

from __future__ import annotations

import csv
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re

BASE = Path(__file__).resolve().parent
WORKSPACE = BASE.parents[1]
OLD = WORKSPACE / "output/bh3_bclw_source_audit_20261020"
CRITERIA = WORKSPACE / "output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json"
EXPECTED_CRITERIA = "23f2be45ab89a718b03079f34f0832d1a6bf8f7789581c6c8c2fd5bad8d22f5b"


def load(path: Path) -> dict:
    """Load an archived UTF-8 JSON object."""
    return json.loads(path.read_text(encoding="utf-8"))


def save(path: Path, obj: object) -> None:
    """Write finite, inspectable research metadata."""
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    """Write an Excel-readable CSV without altering original measurements."""
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_measurement(raw: str) -> dict:
    """Distinguish fit intervals, censored values, and approximate point values."""
    raw = raw.strip()
    match = re.fullmatch(r"([0-9.]+)\s*±\s*([0-9.]+)", raw)
    if match:
        value, halfwidth = map(float, match.groups())
        return {"raw": raw, "status": "finite_fit_CI", "kd_nM": value,
                "reported_95fitCI_halfwidth_nM": halfwidth}
    match = re.fullmatch(r"([<>~])\s*([0-9.]+)[a-z]?", raw)
    if not match:
        raise ValueError(f"Unrecognized measurement {raw!r}")
    symbol, value = match.groups()
    return {"raw": raw, "status": {">": "lower_bound", "<": "upper_bound", "~": "approximate"}[symbol],
            "bound_or_approx_nM": float(value)}


def normalize(name: str) -> str:
    """Keep an explicit alias for the prime-marked clone in Figure 7."""
    return (name.replace("BL6*-16", "BL6prime-16").replace("BL6’-16", "BL6prime-16")
            .replace("_G3eCb", "_G3eC").replace("_I2dRd", "_I2dR"))


def main() -> None:
    """Generate raw measurements, pair audit, sequence graph, and receptor map."""
    if hashlib.sha256(CRITERIA.read_bytes()).hexdigest() != EXPECTED_CRITERIA:
        raise ValueError("Frozen v2 changed")
    article = load(BASE / "extracted/Foight_2015_viral_landscape.json")
    supplement = load(BASE / "extracted/Foight2015_supplement1.json")
    tables = supplement["tables"]
    assert tables[5][0][5] == "Bcl-w" and len(tables[5]) == 30
    assert len(tables[6]) == 25
    raw_records = []
    measurements = {}
    for row in tables[5][1:]:
        name = normalize(row[0])
        measurement = parse_measurement(row[5])
        measurements[name] = measurement
        raw_records.append({"source": "Foight2015", "table": "Supplementary Table 6",
                            "receptor": "Bcl-w", "peptide": name, "unit": "nM", **measurement})
    table1 = article["tables"][0]
    row = next(row for row in table1["rows"] if row and row[0] == "Bcl-w")
    for name, raw in zip(["Bim", "Bim_R3bS", "Bim_A2eT", "Bim_A2eT_E2gG"], row[1:], strict=True):
        measurement = parse_measurement(raw)
        measurements[name] = measurement
        raw_records.append({"source": "Foight2015", "table": "Table 1", "receptor": "Bcl-w",
                            "peptide": name, "unit": "nM", **measurement})
    # A fixed union of columns preserves censored values without point imputation.
    columns = ["source", "table", "receptor", "peptide", "unit", "raw", "status", "kd_nM",
               "reported_95fitCI_halfwidth_nM", "bound_or_approx_nM"]
    write_csv(BASE / "foight2015_bclw_raw_Kd.csv", [{c: r.get(c, "") for c in columns} for r in raw_records])

    clone_sequences = {normalize(row[0]): row[1] for row in tables[3][1:] if len(row) == 3 and row[0]}
    # Alphabetical names are deep-sequencing-selected sequences, not source-designated parents.
    figure7_alias_sequences = {
        "KL_a": "RPHIWMHQGLRRLGDELNAYYAR", "KL_b": "RPKIWIDQGLRRLGDELNAYYAR",
        "BL_a": "RPIIWIAQQISRVGDEFNAYYAR", "BL_b": "RPIIWIAQQIGRVGDRENAYYAR",
        "BL_c": "RPIIWIAQQIGRVGDEENAYYAR", "BL_d": "RPIIWIAQQIGRVGDRFNAYYAR",
    }
    sequences = dict(figure7_alias_sequences)
    original_names = [normalize(row[0]) for row in tables[5][1:23]]
    for name in original_names:
        if name in sequences:
            continue
        library = clone_sequences[name]
        assert len(library) == 25 and library.endswith("RV"), (name, library)
        sequences[name] = library[:-2]
    previous_sequences = {r["name"]: r["sequence"] for r in load(OLD / "peptide_sequences_source.json")["records"]}
    for name in ["Bim", "Bim_A2eT", "Bim_A2eT_E2gG"]:
        sequences[name] = previous_sequences[name]
    sequences["Bim_R3bS"] = sequences["Bim"][:10] + "S" + sequences["Bim"][11:]
    named_mutants = [
        ("KL6-7", "G3eC", 14, "G", "C"), ("KL6-7", "Y4eK", 21, "Y", "K"),
        ("KL6-7", "D3fA", 15, "D", "A"), ("KL6-7", "I2dR", 6, "I", "R"),
        ("BL6-22", "G3eC", 14, "G", "C"), ("BL6-22", "Y4eK", 21, "Y", "K"),
    ]
    for parent, mutation, pos, wild, altered in named_mutants:
        assert sequences[parent][pos - 1] == wild
        sequences[f"{parent}_{mutation}"] = sequences[parent][:pos - 1] + altered + sequences[parent][pos:]
    sequences["KScomp"] = article["tables"][1]["rows"][5][1]
    assert set(sequences) == set(measurements), (set(sequences) - set(measurements), set(measurements) - set(sequences))
    assert all(len(s) == 23 for s in sequences.values()), {n: len(s) for n, s in sequences.items() if len(s) != 23}
    save(BASE / "foight2015_assay_sequences.json", {
        "source": "Foight2015 Figure 7 (PDF page 39 / printed manuscript page 38), Supplementary Table 4, main Tables 1-2, and named mutation definitions",
        "all_peptides": "23 standard amino acids; N-terminal 5/6-fluorescein amidite and C-terminal amide",
        "library_to_assay_note": "Library Table S4 strings are 25 aa ending RV; Figure 7 and Methods specify the 23-aa synthesized assay molecules. Trailing RV is absent in assay sequences.",
        "alias": {"BL6*-16": "BL6prime-16", "BL6’-16": "BL6prime-16",
                  "KL6-7_G3eCb": "KL6-7_G3eC", "BL6-22_G3eCb": "BL6-22_G3eC",
                  "KL6-7_I2dRd": "KL6-7_I2dR"},
        "sequences": sequences,
        "sequence_source_note": "Six alphabetically named sequences were transcribed and visually checked from Figure 7; their names do not establish mutation direction."
    })
    edges = [
        ("Bim", "Bim_R3bS", "Table 1; named R3bS mutant", False),
        ("Bim", "Bim_A2eT", "Table 1; named A2eT mutant", True),
        ("Bim_A2eT", "Bim_A2eT_E2gG", "Table 1; named E2gG mutant", True),
        ("BL6-19", "BL_d", "Main text explicitly describes E3g-to-arginine substitution", False),
        ("BL_c", "BL_b", "Main text explicitly describes E3g-to-arginine substitution", False),
    ]
    edges += [(p, f"{p}_{m}", "Named single mutant in Table S6", False) for p, m, *_ in named_mutants]
    pair_rows = []
    for parent, variant, basis, duplicate in edges:
        differences = [(i + 1, a, b) for i, (a, b) in enumerate(zip(sequences[parent], sequences[variant], strict=True)) if a != b]
        assert len(differences) == 1, (parent, variant, differences)
        pos, wild, altered = differences[0]
        p, v = measurements[parent], measurements[variant]
        ratio, low, high = "", "", ""
        point = "not_identifiable_both_censored"
        if p["status"] == v["status"] == "finite_fit_CI":
            ratio = v["kd_nM"] / p["kd_nM"]
            point = "gain_candidate" if ratio <= .5 else "loss_candidate" if ratio >= 2 else "near_neutral"
            low = (v["kd_nM"] - v["reported_95fitCI_halfwidth_nM"]) / (p["kd_nM"] + p["reported_95fitCI_halfwidth_nM"])
            high = (v["kd_nM"] + v["reported_95fitCI_halfwidth_nM"]) / (p["kd_nM"] - p["reported_95fitCI_halfwidth_nM"])
        pair_rows.append({"parent": parent, "variant": variant, "mutation_local_1based": f"{wild}{pos}{altered}",
                          "orientation_basis": basis, "parent_Kd_nM_reported": p["raw"],
                          "variant_Kd_nM_reported": v["raw"], "ratio_variant_over_parent": ratio,
                          "point_candidate": point, "marginal_endpoint_envelope_low_NOT_joint95": low,
                          "marginal_endpoint_envelope_high_NOT_joint95": high,
                          "edge_identity_already_in_Foight2014": duplicate,
                          "formal_v2_label_released": False})
    assert len(pair_rows) == 11
    write_csv(BASE / "foight2015_source_directed_pair_audit.csv", pair_rows)

    unordered = []
    for first, second in combinations(sorted(original_names), 2):
        differences = [(i + 1, a, b) for i, (a, b) in enumerate(zip(sequences[first], sequences[second], strict=True)) if a != b]
        if len(differences) != 1:
            continue
        pos, a, b = differences[0]
        first_kd, second_kd = measurements[first], measurements[second]
        finite = first_kd["status"] == second_kd["status"] == "finite_fit_CI"
        unordered.append({"first_name_lexical_no_parent_claim": first, "second_name_lexical_no_parent_claim": second,
                          "local_position_1based": pos, "first_aa": a, "second_aa": b,
                          "first_Kd_nM_reported": first_kd["raw"], "second_Kd_nM_reported": second_kd["raw"],
                          "point_fold_difference_if_both_finite": max(first_kd["kd_nM"], second_kd["kd_nM"]) / min(first_kd["kd_nM"], second_kd["kd_nM"]) if finite else "",
                          "source_directed_edge_in_previous_table": any({r["parent"], r["variant"]} == {first, second} for r in pair_rows),
                          "count_as_new_v2_gain": False,
                          "note": "Sequence-neighbour audit only. No parent is inferred from Kd or alphabetic order."})
    write_csv(BASE / "foight2015_unordered_single_substitution_contrasts.csv", unordered)

    canonical = load(OLD / "receptor_core_audit.json")["Bcl-w"]["native_core"]
    published_core = next(r[1] for r in tables[6] if r[0] == "Bcl-w")
    published_full = next(r[1] for r in load(BASE / "extracted/DeBartolo2014_TableS6_constructs.json")["tables"][0] if r[0] == "myc-Bcl-w SPOT array & FP")
    diffs = [{"core_position_1based": i + 1, "UniProt_Q92843_canonical": a, "published_construct": b}
             for i, (a, b) in enumerate(zip(canonical, published_core, strict=True)) if a != b]
    assert diffs == [{"core_position_1based": 133, "UniProt_Q92843_canonical": "Q", "published_construct": "R"}]
    assert published_full.endswith(published_core)
    prefix = published_full[:-len(published_core)]
    assert prefix == "GSEQKLISEEDLGGSGGTS" and len(published_full) == 183
    variants = [f for f in load(BASE / "raw/Q92843.json")["features"] if f.get("featureId") == "VAR_048418"]
    assert len(variants) == 1 and variants[0]["alternativeSequence"]["alternativeSequences"] == ["R"]
    save(BASE / "receptor_construct_reconciliation.json", {
        "status": "PUBLIC_CONSTRUCT_SEQUENCE_RESOLVED_FOR_DEBARTOLO2014_AND_FOIGHT2015_NOT_DIRECTLY_CONFIRMED_FOR_FOIGHT2014_BATCH",
        "published_core_sequence": published_core, "core_length": len(published_core),
        "published_post_TEV_prefix": prefix, "prefix_length": len(prefix),
        "post_TEV_full_sequence": published_full, "post_TEV_full_length": len(published_full),
        "source_sequence_agreement": "DeBartolo2014 Table S6 full myc-Bcl-w matches Foight2015 Table S7 core plus the post-TEV tag inferred from its pSVM tag sequence.",
        "canonical_differences": diffs, "UniProt_variant": variants[0],
        "Foight2014_direct_evidence_limit": "Foight2014 refers to earlier expression methods but does not directly give this clone sequence. Transferring R133 and the 19-aa prefix to its exact assay batch is an inference, not a verified batch identity.",
        "effect_on_existing_input": "The prepared canonical-Q133 core-only FASTA omits both the tag and the R133 sequence in related published constructs. It remains an explicit unverified proxy and must not be called exact assay input. No silent replacement or new prediction was performed.",
        "biological_effect": "No affinity effect of receptor Q133R is asserted or computed. Natural-variant annotation is not proof of unchanged binding."
    })
    existing_rows = list(csv.DictReader((OLD / "candidate_source_audit.csv").open(encoding="utf-8-sig")))
    bclw = [r for r in existing_rows if r["receptor"] == "Bcl-w"]
    new_losses = [r for r in pair_rows if r["point_candidate"] == "loss_candidate" and not r["edge_identity_already_in_Foight2014"]]
    assert len(new_losses) == 2 and sum(r["point_candidate"] == "gain_candidate" for r in pair_rows) == 0
    summary = {
        "date": "2026-10-20", "status": "BOUNDED_SOURCE_AUDIT_COMPLETE_NO_NEW_MODEL_COMPUTATION",
        "frozen_v2_sha256": EXPECTED_CRITERIA,
        "existing_Foight2014_Bclw": {"source_directed_pairs": len(bclw), "gain_point_candidates": 2,
                                    "loss_point_candidates": 5, "source_parents": sorted({r["parent"] for r in bclw}),
                                    "gain_parent_backgrounds": ["Bim_A2eT"], "formal_labels_released": 0},
        "Foight2015_Bclw": {"raw_measurements": len(raw_records), "source_directed_single_substitution_edges": len(pair_rows),
                            "source_directed_gain_point_candidates": 0,
                            "source_directed_loss_point_candidates": sum(r["point_candidate"] == "loss_candidate" for r in pair_rows),
                            "new_edge_identity_loss_candidates": len(new_losses),
                            "both_censored_source_directed_pairs": sum(r["point_candidate"] == "not_identifiable_both_censored" for r in pair_rows),
                            "original_library_peptides": len(original_names), "unordered_single_substitution_contrasts": len(unordered),
                            "formal_labels_released": 0},
        "full_panel_decision": "HOLD_NO_QUALIFYING_BCLW_EXPANSION",
        "reasons": [
            "No new source-directed Bcl-w gain candidate in another background was recovered in this bounded screen.",
            "BL_a and BL6prime-16 are measured one-substitution neighbours (79±11 vs 11±1.3 nM), but the paper selected them from different library sequencing routes without defining parent-to-mutant orientation. Do not choose the higher-Kd sequence as a new parent to manufacture an enhancement.",
            "Overlapping mutation-edge identities across Foight2014 and Foight2015 are not independent new benchmark pairs; assay batches and raw-curve reuse remain unresolved.",
            "Even a purely illustrative union of 7 existing unique edge identities plus 2 new loss candidates gives only 9 pairs and retains just 2 gains concentrated in one parent. This is not an admitted pooled panel.",
            "Public construct sequences reveal R133 versus canonical Q133 and a 19-aa post-TEV prefix; the Foight2014 exact batch mapping remains unconfirmed."
        ],
        "pilot_readiness": "HOLD_EXISTING_INPUT_DRAFT_PENDING_RECEPTOR_IDENTITY_RECONCILIATION_AND_RUN_CONFIGURATION_FREEZE",
        "scope_of_hold": "Bcl-w expansion and current input draft only; not a stop to the peptide-design research project or read-only source preparation.",
        "next_research_action": "End this bounded Bcl-w expansion round; retain the verified raw data and construct map. Prioritize another receptor-assay source with existing source-defined bidirectional Kd point-mutation tables and exact assay inputs. No automatic prediction launch.",
        "new_predictions": 0, "new_scorer_executions": 0, "training_runs": 0,
        "paid_services": 0, "third_party_outreach": 0
    }
    save(BASE / "decision.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
