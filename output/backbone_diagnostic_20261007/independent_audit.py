"""Independent read-only audit of actual structures, positions, tables and statistics.

This script does not import the diagnostic implementation modules.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
BACKBONE = {"N", "CA", "C", "O", "OXT"}
COUNT_TERMS = {"Number of Residues", "Interface Residues", "Interface Residues Clashing", "Interface Residues VdW Clashing", "Interface Residues BB Clashing"}
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))


def js(path):
    return json.loads(path.read_text(encoding="utf-8"))


def rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def close(actual, expected, atol=1e-8):
    assert actual != "", "Missing numeric result"
    actual = float(actual)
    assert math.isfinite(actual) and math.isclose(actual, expected, rel_tol=1e-9, abs_tol=atol), (actual, expected)


def pdb(path, expected):
    atoms, residues = {}, {}
    for line in path.read_text().splitlines():
        assert not line.startswith(("MODEL ", "ENDMDL", "HETATM"))
        if not line.startswith("ATOM  "):
            continue
        assert line[16] == " " and line[26] == " "
        chain, n, atom, residue = line[21], int(line[22:26]), line[12:16].strip(), line[17:20].strip()
        assert chain in expected and residue in AA
        assert (chain, n, atom) not in atoms
        residues[chain, n] = AA[residue]
        xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        assert all(math.isfinite(v) for v in xyz)
        element = line[76:78].strip()
        if element not in ("H", "D") and not atom.startswith(("H", "D")):
            atoms[chain, n, atom] = xyz
    for chain, sequence in expected.items():
        assert sorted(n for c, n in residues if c == chain) == list(range(1, len(sequence) + 1))
        assert "".join(residues[chain, n] for n in range(1, len(sequence) + 1)) == sequence
        assert all((chain, n, atom) in atoms for n in range(1, len(sequence) + 1) for atom in ("N", "CA", "C", "O"))
        assert (chain, len(sequence), "OXT") in atoms
    return atoms


def positions(path):
    value = js(path)
    assert value["units"] == "angstrom" and len(value["atoms"]) == len(value["positions"])
    assert [atom["index"] for atom in value["atoms"]] == list(range(len(value["atoms"])))
    coordinates = {}
    for atom, xyz in zip(value["atoms"], value["positions"]):
        key = atom["chain"], int(atom["residue"]), atom["atom"]
        assert key not in coordinates and len(xyz) == 3 and all(math.isfinite(v) for v in xyz)
        coordinates[key] = tuple(xyz)
    return value, coordinates


def energy_table(path):
    header, records = None, []
    for line in path.read_text().splitlines():
        fields = [v.strip() for v in line.split("\t")]
        if fields[0] == "Pdb":
            assert header is None and len(fields) == len(set(fields))
            header = fields
        elif header is not None and line.strip():
            assert len(fields) == len(header)
            records.append(dict(zip(header, fields)))
    assert len(records) == 1 and header is not None
    r = records[0]
    assert Path(r["Pdb"]).name == "relaxed.pdb" and {r["Group1"], r["Group2"]} == {"A", "B"}
    terms = {k: float(v) for k, v in r.items() if k not in ("Pdb", "Group1", "Group2")}
    assert "Interaction Energy" in terms and all(math.isfinite(v) for v in terms.values())
    return terms


def expected_cases(protocol):
    s = protocol["scientific_protocol"]
    mutation_by_id = {m["peptide_id"]: m for m in s["mutations"]}
    old_pairs = [r for r in rows(ROOT / protocol["source_pair_csv"]) if r["peptide_id"] in mutation_by_id]
    assert len(old_pairs) == 60
    result = {}
    for r in old_pairs:
        assert r["status"] == "ok"
        for role, prefix in (("wt_reference", "wt"), ("mutant", "mutant")):
            uid = f"{r['job_id']}_repeat{r['repeat']}_{role}"
            for condition in s["conditions"]:
                case = dict(r, role=role, source_uid=uid, condition=condition, pair_uid=f"{r['job_id']}_repeat{r['repeat']}",
                            input_sequence=s["parent_sequence"] if role == "wt_reference" else mutation_by_id[r["peptide_id"]]["sequence"],
                            source_file=ROOT / "output/foldx_matched_diagnostic_20261006" / r[f"{prefix}_path"])
                job_id = f"{uid}_{condition}"
                assert job_id not in result
                result[job_id] = case
    assert len(result) == 240
    return result


def direction(gains, expected):
    sign = 1 if expected == "stronger" else -1
    return "correct" if all(v * sign > 0 for v in gains) else "incorrect" if all(v * sign < 0 for v in gains) else "inconclusive"


def unique_table(name, key):
    records = rows(OUT / "analysis" / (name + ".csv"))
    indexed = {tuple(row[field] for field in key): row for row in records}
    assert len(records) == len(indexed), (name, "duplicate row")
    return indexed


def audit_statistics(cases, raw_terms):
    """Reconstruct all pairing and aggregation directly from actual raw tables."""
    structures = unique_table("structure_results", ("job_id",))
    assert set(structures) == {(job_id,) for job_id in cases}
    grouped, comparisons = {}, {}
    for job_id, case in cases.items():
        row = structures[job_id,]
        assert row["status"] == "ok" and row["provenance_valid"] == "True"
        for field in ("source_uid", "pair_uid", "peptide_id", "seed", "model_number", "repeat", "role", "condition"):
            assert row[field] == str(case[field]), (job_id, field)
        assert row["units"] == "kcal/mol"
        close(row["energy_kcal_mol"], case["energy"])
        for field in ("frozen_max_A", "iterations", "iteration_cap"):
            close(row[field], case["relaxation"][field])
        assert row["converged"] == str(case["relaxation"]["converged"])
        grouped.setdefault((case["pair_uid"], case["condition"]), {})[case["role"]] = job_id
    pairs = unique_table("pair_results", ("pair_uid", "condition"))
    assert len(pairs) == 120 and set(pairs) == set(grouped)
    groups = {}
    for key, roles in grouped.items():
        assert set(roles) == {"wt_reference", "mutant"}
        wt, mut = cases[roles["wt_reference"]], cases[roles["mutant"]]
        gain = wt["energy"] - mut["energy"]
        row = pairs[key]
        assert row["status"] == "ok" and row["units"] == "kcal/mol"
        assert row["wt_job_id"] == roles["wt_reference"] and row["mutant_job_id"] == roles["mutant"]
        close(row["wt_energy_kcal_mol"], wt["energy"])
        close(row["mutant_energy_kcal_mol"], mut["energy"])
        close(row["gain_kcal_mol"], gain)
        comparisons.setdefault(key[0], {})[key[1]] = gain
        group = (mut["peptide_id"], mut["condition"], mut["seed"], mut["model_number"])
        groups.setdefault(group, {})[int(mut["repeat"])] = gain
    contrasts = unique_table("contrasts", ("pair_uid",))
    assert len(contrasts) == 60 and set(contrasts) == {(uid,) for uid in comparisons}
    contrast_groups = {}
    for uid, values in comparisons.items():
        assert set(values) == {"fixed", "flexible"}
        row = contrasts[uid,]
        assert row["status"] == "ok"
        difference = values["flexible"] - values["fixed"]
        close(row["fixed_gain_kcal_mol"], values["fixed"])
        close(row["flexible_gain_kcal_mol"], values["flexible"])
        close(row["flexible_minus_fixed_gain_kcal_mol"], difference)
        key = (row["peptide_id"], row["seed"], row["model_number"])
        contrast_groups.setdefault(key, {})[int(row["repeat"])] = difference
    backbones = unique_table("backbone_summary", ("peptide_id", "condition", "seed", "model_number"))
    assert len(backbones) == 24 and set(backbones) == set(groups)
    seed_groups = {}
    for key, repeats in groups.items():
        assert set(repeats) == set(range(1, 6))
        values = list(repeats.values())
        row = backbones[key]
        assert row["complete"] == "True" and row["status"] == "complete"
        for field, truth in dict(n_complete_pairs=5, median_gain_kcal_mol=statistics.median(values),
                                 minimum_gain_kcal_mol=min(values), maximum_gain_kcal_mol=max(values),
                                 gain_range_kcal_mol=max(values) - min(values), positive_pairs=sum(v > 0 for v in values),
                                 negative_pairs=sum(v < 0 for v in values), zero_pairs=sum(v == 0 for v in values)).items():
            close(row[field], truth)
        ckey = (key[0], key[2], key[3])
        assert set(contrast_groups[ckey]) == set(range(1, 6))
        close(row["median_flexible_minus_fixed_gain_kcal_mol"], statistics.median(contrast_groups[ckey].values()))
        seed_groups.setdefault(key[:3], {})[int(key[3])] = statistics.median(values)
    seeds = unique_table("seed_summary", ("peptide_id", "condition", "seed"))
    assert len(seeds) == 8 and set(seeds) == set(seed_groups)
    directions = []
    for key, models in seed_groups.items():
        assert set(models) == {1, 2, 3}
        values, row = list(models.values()), seeds[key]
        assert row["complete"] == "True"
        for field, truth in dict(n_complete_backbones=3, median_gain_kcal_mol=statistics.median(values),
                                 minimum_backbone_gain_kcal_mol=min(values), maximum_backbone_gain_kcal_mol=max(values),
                                 positive_backbones=sum(v > 0 for v in values), negative_backbones=sum(v < 0 for v in values),
                                 zero_backbones=sum(v == 0 for v in values)).items():
            close(row[field], truth)
        expected = "stronger" if key[0] == "P53_P27A" else "weaker"
        assert row["experimental_direction"] == expected
        assert row["diagnostic_direction"] == direction(values, expected)
        close(row["median_flexible_minus_fixed_gain_kcal_mol"], statistics.median(
            statistics.median(contrast_groups[key[0], key[2], str(model)].values()) for model in (1, 2, 3)))
        directions.append({field: row[field] for field in ("peptide_id", "condition", "seed", "diagnostic_direction")})
    terms = unique_table("energy_terms", ("job_id", "table_path", "term"))
    expected_terms = {key: value for key, value in raw_terms.items() if Path(key[1]).name.startswith("Interaction_")}
    assert set(terms) == set(expected_terms)
    for key, value in expected_terms.items():
        close(terms[key]["value"], value)
        assert terms[key]["units"] == ("count" if key[2] in COUNT_TERMS else "kcal/mol")
    for filename, subset, count, pair_count in (("pilot_decision.json", [v for v in cases.values() if int(v["repeat"]) == 1], 48, 24),
                                               ("decision.json", list(cases.values()), 240, 120)):
        decision = js(OUT / "analysis" / filename)
        assert decision["technical_pass"] is True and decision["gate_uses_score_direction"] is False
        for field, truth in dict(expected_structures=count, observed_structures=count, valid_structures=count,
                                 expected_pairs=pair_count, complete_pairs=pair_count,
                                 converged_structures=sum(v["relaxation"]["converged"] for v in subset),
                                 nonconverged_structures=sum(not v["relaxation"]["converged"] for v in subset)).items():
            assert decision[field] == truth, (filename, field)
        assert not decision["failure_reasons"]
    decision = js(OUT / "analysis/decision.json")
    assert decision["complete"] is True and decision["complete_contrasts"] == 60 and decision["biological_replicates"] is False
    assert sorted([dict(r, seed=str(r["seed"])) for r in decision["diagnostic_directions"]],
                  key=lambda r: (r["peptide_id"], r["condition"], r["seed"])) == sorted(
        directions, key=lambda r: (r["peptide_id"], r["condition"], r["seed"]))
    assert js(OUT / "analysis/pilot_decision.json")["expand_to_full"] is True
    return dict(structure_rows=240, condition_pair_rows=120, paired_contrast_rows=60, backbone_rows=24, seed_rows=8,
                published_interaction_terms=len(terms), independently_recomputed_directions=directions)


def audit_baseline(protocol):
    inspected = js(OUT / "existing_structure_inspection.json")
    geometry = unique_table("baseline_geometry", ("peptide_id", "seed", "model_number"))
    torsions = unique_table("baseline_torsions", ("peptide_id", "seed", "model_number", "residue"))
    hbonds = unique_table("baseline_hbond_proxies", ("peptide_id", "seed", "model_number", "acceptor_residue", "donor_residue"))
    assert len(geometry) == 18 and len(torsions) == 216 and len(hbonds) == 144
    for item in inspected["structures"]:
        key = item["peptide_id"], str(item["seed"]), str(item["model"])
        row = geometry[key]
        assert sha(ROOT / row["source_path"]) == row["source_sha256"]
        close(row["alpha_like_8_11"], item["alpha_like_B8_11"])
        for angle in item["torsions"]:
            published = torsions[key + (str(angle["position"]),)]
            assert published["alpha_like"] == str(angle["alpha_like"])
            for field in ("phi_deg", "psi_deg"):
                if angle[field] is None:
                    assert published[field] == "" and published[field.replace("_deg", "_status")] == "terminal_missing"
                else:
                    close(published[field], angle[field], 2e-4)
        close(hbonds[key + ("7", "11")]["distance_A"], item["O7_N11_distance_A"], 2e-5)
        assert hbonds[key + ("7", "11")]["eligible_donor"] == str(item["N11_is_possible_amide_donor"])
    for item in inspected["aligned_comparisons"]:
        row = geometry[item["peptide_id"], str(item["seed"]), str(item["model"])]
        for field in ("A_backbone_rmsd_A", "B_backbone_rmsd_A", "B8_12_backbone_rmsd_A"):
            close(row[field], item[field], 2e-5)
    old = {(r["job_id"], r["repeat"], r["role"], r["term"]): float(r["value"]) for r in rows(ROOT / protocol["source_energy_terms_csv"])
           if Path(r["table_path"]).name.startswith("Interaction_") and r["peptide_id"] in ("P53_P27A", "P53_L22A")}
    published = unique_table("baseline_energy_term_contrasts", ("job_id", "repeat", "term"))
    assert len(published) * 2 == len(old)
    for (job, repeat, term), row in published.items():
        wt, mutant = old[job, repeat, "wt_reference", term], old[job, repeat, "mutant", term]
        close(row["wt_reference"], wt)
        close(row["mutant"], mutant)
        close(row["wt_minus_mutant"], wt - mutant)
        assert row["units"] == ("count" if term in COUNT_TERMS else "kcal/mol")
    return dict(structures=18, torsion_rows=216, distance_proxy_rows=144, original_term_contrasts=len(published))


def audit_job_geometry(initial, final, mapping, sequence, relaxation):
    from Bio.PDB.vectors import Vector, calc_dihedral
    derived = {"torsions": [], "hbond_proxies": [], "residue_displacements": []}
    for moment, coordinates in (("before", initial), ("after", final)):
        recorded_torsions = relaxation["torsions_" + moment]
        recorded_hbonds = relaxation["hbond_proxies_" + moment]
        assert len(recorded_torsions) == 12 and len(recorded_hbonds) == 8
        for n in range(1, 13):
            def angle(keys):
                return math.degrees(calc_dihedral(*(Vector(*coordinates[key]) for key in keys)))
            phi = angle([("B", n - 1, "C"), ("B", n, "N"), ("B", n, "CA"), ("B", n, "C")]) if n > 1 else None
            psi = angle([("B", n, "N"), ("B", n, "CA"), ("B", n, "C"), ("B", n + 1, "N")]) if n < 12 else None
            row = recorded_torsions[n - 1]
            assert row["residue"] == n and row["amino_acid"] == sequence[n - 1]
            for field, value in (("phi_deg", phi), ("psi_deg", psi)):
                if value is None:
                    assert row[field] is None and row[field.replace("_deg", "_status")] == "terminal_missing"
                else:
                    close(row[field], value, 1e-6)
            alpha = phi is not None and psi is not None and -100 <= phi <= -30 and -80 <= psi <= -5
            assert row["alpha_like"] == alpha and row["interior_2_11"] == (2 <= n <= 11) and row["local_8_11"] == (8 <= n <= 11)
            derived["torsions"].append(dict(row, moment=moment))
        for n in range(1, 9):
            distance = math.dist(coordinates["B", n, "O"], coordinates["B", n + 4, "N"])
            row = recorded_hbonds[n - 1]
            assert (row["acceptor_residue"], row["donor_residue"]) == (n, n + 4)
            close(row["distance_A"], distance)
            assert row["eligible_donor"] == (sequence[n + 3] != "P")
            assert row["proxy_present"] == (sequence[n + 3] != "P" and distance <= 3.5)
            derived["hbond_proxies"].append(dict(row, moment=moment))
    recorded = relaxation["residue_displacements"]
    assert len(recorded) == 12
    for n in range(1, 13):
        row = recorded[n - 1]
        assert row["residue"] == n and row["amino_acid"] == sequence[n - 1]
        for prefix, backbone_only in (("backbone", True), ("all_heavy", False)):
            keys = [(a["chain"], a["residue"], a["atom"]) for a in mapping if a["chain"] == "B" and a["residue"] == n
                    and a["element"] not in ("H", "D") and (not backbone_only or a["atom"] in BACKBONE)]
            distances = [math.dist(initial[key], final[key]) for key in keys]
            close(row[prefix + "_rmsd_A"], math.sqrt(statistics.mean(d ** 2 for d in distances)))
            close(row[prefix + "_max_A"], max(distances))
            assert row[prefix + "_n_atoms"] == len(keys)
        derived["residue_displacements"].append(row)
    after_angles = [row for row in derived["torsions"] if row["moment"] == "after"]
    assert relaxation["alpha_like_2_11"] == sum(row["alpha_like"] for row in after_angles if 2 <= row["residue"] <= 11)
    assert relaxation["alpha_like_8_11"] == sum(row["alpha_like"] for row in after_angles if 8 <= row["residue"] <= 11)
    assert relaxation["hbond_distance_proxy_count"] == sum(row["proxy_present"] for row in derived["hbond_proxies"] if row["moment"] == "after")
    return derived


def audit_published_geometry(cases):
    counts = {}
    for name, key in (("torsions", ("job_id", "moment", "residue")),
                      ("hbond_proxies", ("job_id", "moment", "acceptor_residue", "donor_residue")),
                      ("residue_displacements", ("job_id", "residue"))):
        published = unique_table(name, key)
        expected = {}
        for job_id, case in cases.items():
            for row in case["geometry"][name]:
                expected[(job_id,) + tuple(str(row[field]) for field in key[1:])] = row
        assert set(published) == set(expected)
        for index, row in expected.items():
            for field, value in row.items():
                if value is None:
                    assert published[index][field] == ""
                elif isinstance(value, bool) or isinstance(value, str):
                    assert published[index][field] == str(value)
                else:
                    close(published[index][field], value)
        counts[name + "_rows"] = len(published)
    return counts


def audit():
    protocol = js(OUT / "protocol.json")
    for name, digest in protocol["source_hashes"].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in protocol["environment"]["external_hashes"].items():
        assert sha(Path(name)) == digest, name
    expected = expected_cases(protocol)
    manifest = rows(OUT / "job_manifest.csv")
    assert len(manifest) == 240 and {r["job_id"] for r in manifest} == set(expected)
    preparation = js(OUT / "preparation.json")
    assert preparation["protocol_sha256"] == sha(OUT / "protocol.json")
    for name, digest in preparation["prepared_hashes"].items():
        assert sha(OUT / name) == digest, name
    for name, digest in preparation["implementation_hashes"].items():
        assert sha(ROOT / name) == digest, name
    state = js(OUT / "execution_state.json")
    assert state["preparation_sha256"] == sha(OUT / "preparation.json")
    assert len(state["jobs"]) == 240 and set(state["jobs"]) == set(expected)
    cases, hydration_keys, raw_terms = {}, {}, {}
    frozen_max, flexible_max, physical_increases = 0., 0., 0
    s = protocol["scientific_protocol"]
    for job_id, case in expected.items():
        job = OUT / "jobs" / job_id
        result = js(job / "result.json")
        tracked = state["jobs"][job_id]
        assert sha(job / "result.json") == (tracked if isinstance(tracked, str) else tracked["result_sha256"])
        assert result["job_id"] == job_id and result["source_uid"] == case["source_uid"]
        assert result["status"] == "ok" and result["provenance_valid"] is True
        actual_artifacts = {p.relative_to(job).as_posix(): sha(p) for p in job.rglob("*")
                            if p.is_file() and p.name != "result.json" and "attempts" not in p.relative_to(job).parts}
        assert actual_artifacts == result["artifact_hashes"], (job_id, "artifact completeness")
        for name, digest in result["artifact_hashes"].items():
            assert sha(job / name) == digest, (job_id, name)
        expected_sequence = {"A": s["receptor_sequence"], "B": case["input_sequence"]}
        source = pdb(case["source_file"], expected_sequence)
        common = OUT / "sources" / case["source_uid"] / "common"
        initial_json, initial = positions(common / "positions.json")
        final_json, final = positions(job / "positions.json")
        assert initial_json["atoms"] == final_json["atoms"] and set(initial) == set(final)
        previous = hydration_keys.setdefault(case["source_uid"], sha(common / "positions.json"))
        assert previous == sha(common / "positions.json")
        hydration = js(common / "hydration.json")
        assert hydration["source_uid"] == case["source_uid"]
        for name, digest in hydration["artifact_hashes"].items():
            assert sha(common / name) == digest, (job_id, "common", name)
        for key, xyz in source.items():
            assert key in initial and math.dist(initial[key], xyz) <= 1e-6, (job_id, key, "hydration moved heavy atom")
        relaxed = pdb(job / "relaxed.pdb", expected_sequence)
        assert set(relaxed) == set(source)
        for key, xyz in relaxed.items():
            assert max(abs(a - b) for a, b in zip(xyz, final[key])) <= .000501, (job_id, key, "PDB serialization")
        frozen = [index for index, atom in enumerate(initial_json["atoms"]) if atom["atom"] in BACKBONE and (atom["chain"] == "A" or case["condition"] == "fixed")]
        rel = js(job / "relaxation.json")
        assert rel["condition"] == case["condition"] and rel["frozen_atom_indices"] == frozen
        assert result["relaxation"] == rel
        assert rel["hydration_sha256"] == sha(common / "hydration.json")
        assert rel["start_positions_sha256"] == sha(common / "positions.json")
        assert final_json["start_positions_sha256"] == sha(common / "positions.json")
        assert rel["platform"] == "CPU" and rel["platform_properties"] == {"Threads": "1", "DeterministicForces": "true"}
        assert rel["iteration_cap"] == 2000 and 0 <= rel["iterations"] <= 2000
        assert rel["iterations"] == len(rel["reporter"])
        assert rel["converged"] == (rel["after"]["free_force_rms_kJ_mol_nm"] <= 10)
        for moment in ("before", "after"):
            assert all(math.isfinite(v) for v in rel[moment].values())
            close(rel[moment]["free_force_rms_kJ_mol_nm"], rel[moment]["free_force_component_rms_kJ_mol_nm"] * math.sqrt(3), 1e-7)
        before_total = sum(rel["before"][name] for name in ("physical_potential_kJ_mol", "tether_potential_kJ_mol"))
        after_total = sum(rel["after"][name] for name in ("physical_potential_kJ_mol", "tether_potential_kJ_mol"))
        assert after_total <= before_total + max(1e-4, abs(before_total) * 1e-9), (job_id, "minimization objective increased")
        physical_increases += rel["after"]["physical_potential_kJ_mol"] > rel["before"]["physical_potential_kJ_mol"]
        frozen_distances = [math.dist(initial[key], final[key]) for key in initial if key[2] in BACKBONE and (key[0] == "A" or case["condition"] == "fixed")]
        largest = max(frozen_distances)
        assert largest <= s["frozen_atom_tolerance_A"], (job_id, largest)
        close(rel["frozen_max_A"], largest)
        frozen_max = max(frozen_max, largest)
        bdist = [math.dist(initial[key], final[key]) for key in initial if key[0] == "B" and key[2] in BACKBONE]
        expected_tether_kJ_mol = .5 * .1 * 4.184 * sum(distance ** 2 for distance in bdist) if case["condition"] == "flexible" else 0.
        close(rel["before"]["tether_potential_kJ_mol"], 0., 1e-8)
        close(rel["after"]["tether_potential_kJ_mol"], expected_tether_kJ_mol, 1e-5)
        close(rel["B_backbone_max_A"], max(bdist))
        close(rel["B_backbone_rmsd_A"], math.sqrt(statistics.mean(distance ** 2 for distance in bdist)))
        if case["condition"] == "flexible":
            flexible_max = max(flexible_max, max(bdist))
        assert js(job / "analyse.command.json")["returncode"] == 0
        interaction = list(job.glob("Interaction_*.fxout"))
        summary = list(job.glob("Summary_*.fxout"))
        assert len(interaction) == len(summary) == 1
        it, st = energy_table(interaction[0]), energy_table(summary[0])
        assert "IntraclashesGroup1" in it and "IntraclashesGroup2" in it
        close(it["Interaction Energy"], st["Interaction Energy"], 1e-6)
        close(result["raw"]["energy_kcal_mol"], it["Interaction Energy"])
        assert result["raw"]["tables"] == {interaction[0].name: it, summary[0].name: st}
        for path, terms in ((interaction[0], it), (summary[0], st)):
            for term, value in terms.items():
                raw_terms[job_id, path.relative_to(OUT).as_posix(), term] = value
        geometry = audit_job_geometry(initial, final, initial_json["atoms"], case["input_sequence"], rel)
        cases[job_id] = dict(case, energy=it["Interaction Energy"], peptide_max_A=max(bdist), relaxation=rel, geometry=geometry)
    assert len(hydration_keys) == 120 and flexible_max > .01
    assert set(state["hydrations"]) == set(hydration_keys) and not state["invalid_sources"]
    for uid, digest in state["hydrations"].items():
        assert sha(OUT / "sources" / uid / "common/hydration.json") == digest
    published_hashes = js(OUT / "analysis/published_artifact_hashes.json")
    for name, digest in published_hashes.items():
        assert sha(OUT / name) == digest, name
    statistics_audit = audit_statistics(cases, raw_terms)
    baseline_audit = audit_baseline(protocol)
    geometry_audit = audit_published_geometry(cases)
    return dict(status="passed", frozen_workspace_sources=len(protocol["source_hashes"]), frozen_external_sources=len(protocol["environment"]["external_hashes"]),
                structures=len(cases), common_hydrated_inputs=len(hydration_keys), raw_tables=480, raw_numeric_terms=len(raw_terms),
                maximum_frozen_atom_displacement_A=frozen_max, maximum_flexible_peptide_displacement_A=flexible_max,
                physical_potential_increased_structures=physical_increases, statistics=statistics_audit, baseline=baseline_audit, geometry=geometry_audit,
                audit_script_sha256=sha(Path(__file__)), limitation="Correct computation is not affinity/free-energy validation"), cases, raw_terms


if __name__ == "__main__":
    result, _, _ = audit()
    (OUT / "analysis/independent_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
