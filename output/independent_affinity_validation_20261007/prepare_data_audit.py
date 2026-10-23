"""Prepare literature/input audits only. No network, predictions or scores.

Run with the existing graphpep_env Python in WSL. All writes stay in this
new output directory; raw sources and historical results are read-only.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from html.parser import HTMLParser
from pathlib import Path

from Bio.PDB.MMCIF2Dict import MMCIF2Dict

TASK = Path(__file__).resolve().parent
ROOT = TASK.parent.parent
SOURCES = TASK / "sources"
AA = set("ACDEFGHIKLMNPQRSTVWY")
PAPER = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8463443/"
BASE = "TSFAEYWNLASP"
PANELS = {
    "MDM2": [("L10A", BASE), ("T1L", "LSFAEYWNLASP"),
             ("A4L", "TSFLEYWNLASP"), ("N8A", "TSFAEYWALASP"),
             ("P12Q", "TSFAEYWNLASQ")],
    "MDMX": [("L10A", BASE), ("T1L", "LSFAEYWNLASP"),
             ("A4L", "TSFLEYWNLASP"), ("N8A", "TSFAEYWALASP"),
             ("S11M", "TSFAEYWNLAMP")],
}
PRIMARY_EDGES = {
    "MDM2": [("L10A", "A4L", "stronger", "forward"),
             ("L10A", "N8A", "stronger", "forward"),
             ("T1L", "L10A", "weaker", "reverse_of_published_enhancement"),
             ("P12Q", "L10A", "weaker", "reverse_of_published_enhancement")],
    "MDMX": [("L10A", "A4L", "stronger", "forward"),
             ("L10A", "N8A", "stronger", "forward"),
             ("T1L", "L10A", "weaker", "reverse_of_published_enhancement"),
             ("S11M", "L10A", "weaker", "reverse_of_published_enhancement")],
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(name: str, obj) -> None:
    (TASK / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(name: str, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Empty output: {name}")
    with (TASK / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class TableReader(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tables = []
        self.current = self.row = self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self.current = []
        elif self.current is not None and tag == "tr":
            self.row = []
        elif self.row is not None and tag in ("td", "th"):
            self.cell = []

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self.cell is not None:
            self.row.append("".join(self.cell).strip())
            self.cell = None
        elif tag == "tr" and self.row is not None:
            self.current.append(self.row)
            self.row = None
        elif tag == "table" and self.current is not None:
            self.tables.append(self.current)
            self.current = None


def parse_value(text: str) -> tuple[float, float]:
    matched = re.fullmatch(r"\s*([\d.]+)\s*±\s*([\d.]+)\s*", text)
    if matched is None:
        raise ValueError(f"Unexpected reported quantity: {text!r}")
    return tuple(map(float, matched.groups()))


def project_exposure(sequences: set[str]) -> tuple[dict, list]:
    found = {sequence: [] for sequence in sequences}
    inspected = []
    # Search all historical CSV sequence columns; do not print score values.
    for path in sorted((ROOT / "output").rglob("*.csv")):
        if TASK in path.parents:
            continue
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            columns = [c for c in (reader.fieldnames or []) if c.lower() in
                       {"sequence", "peptide_sequence", "parent_sequence", "child_sequence"}]
            if not columns:
                continue
            inspected.append({"path": str(path.relative_to(ROOT)), "sha256": sha(path), "columns": columns})
            for number, row in enumerate(reader, 2):
                for column in columns:
                    sequence = row.get(column, "").strip()
                    if sequence in found:
                        found[sequence].append({"path": str(path.relative_to(ROOT)),
                                                "line": number, "column": column,
                                                "peptide_id": row.get("peptide_id", "")})
    return found, inspected


def cif_audit(path: Path) -> dict:
    values = MMCIF2Dict(str(path))
    sequences = [re.sub(r"\s+", "", x) for x in values["_entity_poly.pdbx_seq_one_letter_code_can"]]
    entities = dict(zip(values["_entity_poly.entity_id"], sequences))
    chains = []
    auth = values.get("_atom_site.auth_asym_id", [])
    entity_ids = values.get("_atom_site.label_entity_id", [])
    label_ids = values.get("_atom_site.label_seq_id", [])
    atom_names = values.get("_atom_site.label_atom_id", [])
    for chain in sorted(set(auth)):
        indices = [i for i, x in enumerate(auth) if x == chain and entity_ids[i] in entities
                   and label_ids[i].isdigit()]
        if not indices:
            continue
        entity = entity_ids[indices[0]]
        sequence = entities[entity]
        resolved = sorted({int(label_ids[i]) for i in indices})
        missing = sorted(set(range(1, len(sequence) + 1)) - set(resolved))
        chains.append({"auth_chain": chain, "entity_id": entity, "full_sequence": sequence,
                       "full_length": len(sequence), "resolved_residue_count": len(resolved),
                       "missing_label_positions": missing,
                       "resolved_sequence": "".join(sequence[i - 1] for i in resolved),
                       "OXT_label_positions": sorted({int(label_ids[i]) for i in indices if atom_names[i] == "OXT"})})
    dif_keys = ["_struct_ref_seq_dif." + k for k in
                ("pdbx_pdb_strand_id", "seq_num", "mon_id", "db_mon_id", "db_seq_num", "details")]
    differences = [dict(zip([k.split(".")[-1] for k in dif_keys], row))
                   for row in zip(*(values.get(k, []) for k in dif_keys))]
    return {"pdb_id": path.stem, "sha256": sha(path),
            "entity_sequences": entities, "chains": chains,
            "deposited_sequence_differences": differences,
            "resolution_A": values.get("_refine.ls_d_res_high", []),
            "use": "source/construct audit only; no new native scores; not exact structures for primary candidates"}


def local_graphpep_overlap(pdb_ids: list[str]) -> dict:
    graphpep = Path("/home/g/GraphPep_v1.3")
    paths = [graphpep / "Training_set/list", graphpep / "Training_set/split_data/train.csv",
             graphpep / "Training_set/split_data/valid.csv"]
    paths += sorted(graphpep.glob("*_decoys/list*"))
    inspected = []
    matches = {pdb_id.upper(): [] for pdb_id in pdb_ids}
    training_index = {}
    for path in paths:
        if not path.is_file():
            continue
        inspected.append({"path": str(path), "sha256": sha(path), "size_bytes": path.stat().st_size})
        counts = {}
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            token = line.split(",", 1)[0].strip()
            if path.name == "list" and "Training_set" in str(path):
                training_index[token.upper()] = number
            for pdb_id in matches:
                if token[:4].upper() == pdb_id:
                    record = counts.setdefault(token, {"path": str(path), "target": token,
                                                       "first_line": number, "rows": 0})
                    record["rows"] += 1
        for target, record in counts.items():
            matches[target[:4].upper()].append(record)
    write_csv("graphpep_training_targets.csv", [
        {"target_id": target, "list_line": line} for target, line in training_index.items()])
    core = [{"path": str(p), "sha256": sha(p)} for p in sorted((graphpep / "bin").rglob("*.py"))]
    return {"repository": str(graphpep), "inspected_lists": inspected, "pdb_id_matches": matches,
            "training_targets_in_local_list": len(training_index), "core_source_hashes": core,
            "scope": "Exact PDB ID audit of installed lists; no exhaustive homolog/sequence or AFM/ESM/PepMLM training audit",
            "training_clean_claim_allowed": False}


def prepare() -> None:
    if (TASK / "prepared_snapshot.json").exists():
        raise RuntimeError("Existing prepared snapshot must be preserved; prepare a new version instead.")
    criteria = json.loads((TASK / "frozen_criteria_v1.json").read_text(encoding="utf-8"))
    assert criteria["new_scores_observed"] is False
    parser = TableReader()
    parser.feed((SOURCES / "PMC8463443.html").read_text(encoding="utf-8"))
    tables = [t for t in parser.tables if any(len(r) == 9 and r[0] == BASE for r in t)]
    assert len(tables) == 1
    raw = {r[0]: r for r in tables[0] if len(r) == 9 and len(r[0]) == 12 and set(r[0]) <= AA}
    assert len(raw) == 23  # base plus 22 FP-selected analogs; not the 94-member initial library
    literal_rows = []
    for sequence, row in raw.items():
        literal_rows.append(dict(zip(
            ["sequence", "mdm2_kd_reported_nM", "mdm2_inverse_kd_ratio_reported",
             "mdm2_ki_reported_nM", "mdm2_inverse_ki_ratio_reported",
             "mdmx_kd_reported_nM", "mdmx_inverse_kd_ratio_reported",
             "mdmx_ki_reported_nM", "mdmx_inverse_ki_ratio_reported"], row)))
    write_csv("literature_table2.csv", literal_rows)

    all_sequences = {sequence for panel in PANELS.values() for _, sequence in panel}
    exposure, inspected = project_exposure(all_sequences)
    candidates = []
    for target, panel in PANELS.items():
        column = 1 if target == "MDM2" else 5
        for label, sequence in panel:
            kd, error = parse_value(raw[sequence][column])
            candidates.append({"peptide_id": target + "_" + label, "system_id": target + "_PMI_L10A_2021",
                               "target": target, "family": "PMI_L10A_shared_scaffold",
                               "sequence": sequence, "length": len(sequence),
                               "background_reference_id": target + "_L10A",
                               "mutation_from_background": "WT_background" if label == "L10A" else label,
                               "kd_nM": kd, "reported_error_nM": error,
                               "error_type": "not_explicitly_defined_for_SPR_Table2_in_accessed_main_text",
                               "sd_nM": "", "sem_nM": "", "censoring": "none",
                               "assay": "SPR competition", "assay_temperature_C": 25,
                               "assay_buffer": "10 mM HEPES;150 mM NaCl;0.005% P20;pH7.4",
                               "independent_experiments": 3, "technical_duplicates_per_experiment": 2,
                               "receptor_construct_reported": "MDM2(25-109) synthetic" if target == "MDM2" else "MDMX(24-108) synthetic; exact assay WT/AAA unresolved",
                               "N_terminal": "unresolved; free amine is a hypothesis from synthetic lineage",
                               "C_terminal": "unresolved; free acid is a hypothesis from synthetic lineage",
                               "source_doi": "10.1016/j.apsb.2021.06.010", "source_url": PAPER,
                               "source_table": "Table 2; Kd column only; exclude Ki columns",
                               "source_snapshot": "sources/PMC8463443.html",
                               "project_exact_sequence_exposure_count": len(exposure[sequence]),
                               "used_for_project_development": bool(exposure[sequence]),
                               "shared_development_family": True,
                               "exact_candidate_PDB": "not_identified",
                               "related_family_PDB": "5UMM (M3; different peptide)" if target == "MDM2" else "5UML/3EQY (AAA receptor; different peptide)",
                               "prediction_release": "HOLD_data_audit"})
    write_csv("candidate_affinity.csv", candidates)
    index = {row["peptide_id"]: row for row in candidates}
    edges = []
    for target, definitions in PRIMARY_EDGES.items():
        for parent, variant, direction, origin in definitions:
            p, v = index[target + "_" + parent], index[target + "_" + variant]
            differences = [(i + 1, a, b) for i, (a, b) in enumerate(zip(p["sequence"], v["sequence"])) if a != b]
            assert len(differences) == 1
            position, aa_p, aa_v = differences[0]
            ratio = v["kd_nM"] / p["kd_nM"]
            low = (v["kd_nM"] - 2 * v["reported_error_nM"]) / (p["kd_nM"] + 2 * p["reported_error_nM"])
            high = (v["kd_nM"] + 2 * v["reported_error_nM"]) / (p["kd_nM"] - 2 * p["reported_error_nM"])
            strong = high <= 1 / 3 if direction == "stronger" else low >= 3
            assert strong
            edges.append({"pair_id": f"{target}_{parent}_to_{variant}", "system_id": p["system_id"],
                          "parent_id": p["peptide_id"], "variant_id": v["peptide_id"],
                          "parent_sequence": p["sequence"], "variant_sequence": v["sequence"],
                          "peptide_position_1based": position, "mutation": f"{aa_p}{position}{aa_v}",
                          "mutation_count": 1, "expected_direction": direction,
                          "contrast_origin": origin, "kd_variant_over_parent": ratio,
                          "experimental_log10_gain": math.log10(1 / ratio),
                          "reported_error_envelope_ratio_low": low,
                          "reported_error_envelope_ratio_high": high,
                          "error_envelope_is_confidence_interval": False,
                          "threefold_direction_retained_under_reported_error_envelope": strong,
                          "primary_candidate": True, "count_reverse_edge_again": False,
                          "selection_note": "Published panel enriched for FP improvements; reverse contrasts are NOT new loss-of-function assay samples"})
    write_csv("mutation_pairs.csv", edges)

    old_reference = ROOT / "output/affinity_diagnostic_mdm2_20261006/reference_affinity.csv"
    with old_reference.open(encoding="utf-8-sig", newline="") as stream:
        development = [{**row, "dataset_role": "DEVELOPMENT_USED_FOR_METHOD_SELECTION",
                        "eligible_as_new_primary_holdout": False} for row in csv.DictReader(stream)]
    assert len(development) == 10
    write_csv("development_set.csv", development)

    structures = [cif_audit(p) for p in sorted(SOURCES.glob("*.cif"))]
    by_pdb = {s["pdb_id"]: s for s in structures}
    uniprot = {}
    for accession in ("Q00987", "O15151"):
        uniprot[accession] = "".join(line.strip() for line in (SOURCES / (accession + ".fasta")).read_text().splitlines() if not line.startswith(">"))
    mdm2 = uniprot["Q00987"][24:109]
    mdmx = uniprot["O15151"][23:108]
    assert len(mdm2) == len(mdmx) == 85
    assert by_pdb["5UMM"]["entity_sequences"]["1"] == mdm2
    mdmx_deposited = by_pdb["5UML"]["entity_sequences"]["1"]
    mismatches = [{"uniprot_position": i + 24, "WT": a, "deposited": b}
                  for i, (a, b) in enumerate(zip(mdmx, mdmx_deposited)) if a != b]
    assert mismatches == [{"uniprot_position": 68, "WT": "Q", "deposited": "A"},
                           {"uniprot_position": 69, "WT": "Q", "deposited": "A"},
                           {"uniprot_position": 70, "WT": "E", "deposited": "A"}]
    receptor_rows = [
        {"target": "MDM2", "construct": "25-109", "length": 85, "sequence": mdm2,
         "sequence_source": "UniProt Q00987 + exact agreement with 5UMM entity 1 + existing project receptor",
         "assay_sequence_status": "range/protein_identity_consistent; individual reagent mass/terminal confirmation incomplete"},
        {"target": "MDMX", "construct": "24-108_WT", "length": 85, "sequence": mdmx,
         "sequence_source": "UniProt O15151 canonical subsequence",
         "assay_sequence_status": "HOLD_WT_vs_AAA_unresolved"},
        {"target": "MDMX", "construct": "24-108_Q68A_Q69A_E70A", "length": 85, "sequence": mdmx_deposited,
         "sequence_source": "5UML and 3EQY entity 1",
         "assay_sequence_status": "HOLD_do_not_assume_structure_mutations_apply_to_SPR_reagent"}]
    write_csv("receptor_constructs_pending.csv", receptor_rows)
    # Candidate peptide FASTA only, deliberately no executable complex FASTA until chemistry is resolved.
    (TASK / "candidate_peptides_pending.fasta").write_text("".join(
        f">{label}|shared_scaffold|chemistry_UNRESOLVED|DO_NOT_PREDICT\n{sequence}\n"
        for label, sequence in sorted({(label, sequence) for panel in PANELS.values() for label, sequence in panel})), encoding="utf-8")
    jobs = [{"job_id": f"{c['peptide_id']}_seed{seed}_model{model}", "system_id": c["system_id"],
             "peptide_id": c["peptide_id"], "folding_seed": seed, "parameter_model": model,
             "state": "PLANNED_NOT_STARTED_HOLD_DATA_AUDIT"}
            for c in candidates for seed in criteria["folding"]["seeds"] for model in criteria["folding"]["models"]]
    write_csv("planned_jobs.csv", jobs)
    overlap = local_graphpep_overlap(["5UMM", "5UML", "3EQY", "3EQS", "3LNZ", "1YCR", "2PQK", "3KZ0", "3KJ0", "3KJ1", "3KJ2", "3KZE", "4GVD"])
    write_json("graphpep_overlap_audit.json", overlap)
    write_json("structure_source_audit.json", {"structures": structures, "MDMX_deposited_vs_WT": mismatches})
    assert not any(exposure.values()), "Primary candidates already exposed in historical CSVs; stop and revise scope before freezing."
    audit = {
        "audit_kind": "Preparation/source/input audit; NOT new-score or biological verification",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_table2_rows": len(raw), "selected_target_peptide_rows": len(candidates),
        "unique_selected_peptide_sequences": len(all_sequences), "primary_pair_rows": len(edges),
        "planned_structures_total_max": len(jobs), "actually_predicted_structures": 0,
        "actually_scored_new_structures": 0,
        "historical_sequence_csvs_inspected": inspected,
        "exact_candidate_project_exposure": exposure,
        "exposure_scope_limit": "CSV sequence columns in Windows workspace output; not all external WSL history, communications or every model-training corpus",
        "source_numeric_extraction_pass": True, "sequence_length_and_single_mutation_pass": True,
        "unit_and_assay_separation_pass": True, "error_definition_resolved": False,
        "peptide_terminal_chemistry_resolved": False,
        "MDM2_receptor_sequence_consistent": True, "MDMX_exact_assay_receptor_resolved": False,
        "graphpep_training_overlap_known": bool(overlap["pdb_id_matches"]["5UMM"]),
        "GraphPep_AF_or_ESM_training_clean": "not_established",
        "scientific_entry_gate_pass": False, "new_prediction_warranted_now": False,
        "interpretation": "Concrete small-panel preparation is complete; independent-affinity objective validity is not established."
    }
    write_json("input_audit.json", audit)
    old_paths = [ROOT / "PROJECT_HANDOFF_2026-10-07.md",
                 ROOT / "output/backbone_diagnostic_20261007/analysis/interpretation_zh.md",
                 ROOT / "output/foldx_affinity_mdm2_20261006/analysis/report_zh.md",
                 ROOT / "output/affinity_diagnostic_mdm2_20261006/report_zh.md",
                 ROOT / "output/feedback_control_1ycr_20261006/report_zh.md",
                 old_reference, ROOT / "pipeline/affinity_diagnostic.py", ROOT / "pipeline/foldx_affinity.py",
                 ROOT / "pipeline/relax_feedback_panel.py", ROOT / "pipeline/step3_graphpep.py"]
    protected = [{"path": str(p.relative_to(ROOT)), "sha256": sha(p)} for p in old_paths]
    write_json("prepared_snapshot.json", {
        "state": "score_free_preparation_snapshot",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "new_scores_observed": False,
        "criteria_sha256": sha(TASK / "frozen_criteria_v1.json"),
        "prepared_files": {p.name: sha(p) for p in sorted(TASK.iterdir()) if p.is_file() and p.name not in {"prepared_snapshot.json", "research_workplan.md"}},
        "source_files": {p.name: sha(p) for p in sorted(SOURCES.iterdir()) if p.is_file()},
        "historical_source_hashes": protected,
        "runtime": {"preparation_python": sys.version, "GraphPep_repository": "/home/g/GraphPep_v1.3",
                    "graphpep_python": "/home/g/miniconda3/envs/graphpep_env/bin/python",
                    "colabfold_python": "/home/g/miniconda3/envs/colabfold_exp0/bin/python",
                    "FoldX_executable": "/mnt/e/Users/graphpep与pepmlm/tools/foldx/foldx5_Linux/foldx_20261231",
                    "FoldX_executable_sha256": sha(ROOT / "tools/foldx/foldx5_Linux/foldx_20261231"),
                    "model_weights_or_full_env_audit": "not repeated; must freeze before actual prediction"}
    })
    print(json.dumps({"candidate_rows": len(candidates), "pairs": len(edges), "planned_max_structures": len(jobs),
                      "new_predictions": 0, "data_entry_pass": False,
                      "overlap_hits": {k: v for k, v in overlap["pdb_id_matches"].items() if v}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    prepare()
