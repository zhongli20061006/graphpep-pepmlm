"""Frozen, previously viewed DEVELOPMENT diagnostic using matched FoldX controls.

No BuildModel total energy is used as binding energy. Each gain is the
AnalyseComplex Interaction Energy of a specific WT reference minus its mutant.
"""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output/foldx_matched_diagnostic_20261006"
RECEPTOR = "ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVV"
PARENTS = {"P53_WT": "ETFSDLWKLLPE", "PMI_WT": "TSFAEYWNLLSP"}
MUTATIONS = [
    dict(peptide_id="P53_P27A", parent_id="P53_WT", position=11, **{"from": "P", "to": "A"}, sequence="ETFSDLWKLLAE", individual_list="PB11A;", experimental_direction="stronger"),
    dict(peptide_id="P53_L22A", parent_id="P53_WT", position=6, **{"from": "L", "to": "A"}, sequence="ETFSDAWKLLPE", individual_list="LB6A;", experimental_direction="weaker"),
    dict(peptide_id="PMI_N8A", parent_id="PMI_WT", position=8, **{"from": "N", "to": "A"}, sequence="TSFAEYWALLSP", individual_list="NB8A;", experimental_direction="stronger"),
    dict(peptide_id="PMI_W7A", parent_id="PMI_WT", position=7, **{"from": "W", "to": "A"}, sequence="TSFAEYANLLSP", individual_list="WB7A;", experimental_direction="weaker"),
]
SEEDS, MODELS = (20261011, 20261012), (1, 2, 3)
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
SETTINGS = "temperature=298\npH=7.4\nionStrength=0.15\nwater=-IGNORE\nscreen=false\nnoNterm=false\nnoCterm=false\n"
BUILD_CONFIG = "command=BuildModel\npdb=input.pdb\nmutant-file=individual_list.txt\nnumberOfRuns=5\nmoveNeighbours=true\nout-pdb=true\n" + SETTINGS
ANALYSE_TEMPLATE = "command=AnalyseComplex\npdb={pdb}\nanalyseComplexChains=A,B\n" + SETTINGS
PREPARED_NAMES = {"input.pdb", "individual_list.txt", "build.cfg", "analyse_template.cfg"}
MANIFEST_FIELDS = "job_id peptide_id parent_id sequence parent_sequence individual_list experimental_direction seed model_number source_path source_sha256 input_path input_sha256 mutation_sha256 build_config_sha256 analyse_template_sha256".split()
PAIR_FIELDS = "job_id peptide_id parent_id seed model_number repeat experimental_direction status failure_reason provenance_valid job_failure_reason source_path wt_path mutant_path wt_energy_kcal_mol mutant_energy_kcal_mol gain_kcal_mol old_parent_energy_kcal_mol old_mutant_energy_kcal_mol old_independent_gain_kcal_mol wt_repacking_delta_kcal_mol mutant_vs_old_mutant_delta_kcal_mol wt_source_backbone_rmsd_A mutant_source_backbone_rmsd_A pair_backbone_rmsd_A wt_source_backbone_max_A mutant_source_backbone_max_A pair_backbone_max_A backbone_atom_count wt_raw_paths mutant_raw_paths units".split()
TERM_FIELDS = "job_id peptide_id parent_id seed model_number repeat role table_path term value units".split()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def json_text(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n"


def save_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json_text(value), encoding="utf-8")
    temporary.replace(path)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def csv_text(rows: list[dict], fields: list[str] | None = None) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def check_hashes(hashes: dict, root: Path) -> None:
    for name, expected in hashes.items():
        path = root / name
        require(path.is_file() and sha(path) == expected, f"Frozen file missing/changed: {path}")


def implementation_hashes() -> dict:
    return {name: sha(ROOT / name) for name in ("pipeline/foldx_matched_diagnostic.py", "tests/test_foldx_matched_diagnostic.py")}


def validate_pdb(text: str, peptide: str) -> dict:
    """Require numbered complete A/B chains, ordinary residues and free C termini."""
    residues, atoms = {}, {}
    for line in text.splitlines():
        require(not line.startswith(("MODEL ", "ENDMDL")), "Model-delimited PDB unsupported")
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        require(line.startswith("ATOM  ") and len(line) >= 54, "Nonstandard/truncated atom")
        require(line[16] == " " and line[26] == " ", "Altloc/insertion code")
        chain, number, residue, name = line[21], int(line[22:26]), line[17:20].strip(), line[12:16].strip()
        require(chain in ("A", "B") and residue in AA, "Wrong chain/nonstandard residue")
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        element = element or re.sub(r"^[0-9]+", "", name).upper()[:1]
        if element in ("H", "D"):
            continue
        require(element in ("C", "N", "O", "S"), "Unsupported heavy element")
        key = (chain, number)
        require(key not in residues or residues[key] == residue, "Conflicting residue identity")
        residues[key] = residue
        key = (chain, number, name)
        require(key not in atoms, "Duplicate atom")
        coords = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        require(all(math.isfinite(c) for c in coords), "Nonfinite coordinates")
        atoms[key] = coords
    for chain, sequence in (("A", RECEPTOR), ("B", peptide)):
        numbers = sorted(n for c, n in residues if c == chain)
        require(numbers == list(range(1, len(sequence) + 1)), f"Incomplete/noncontiguous chain {chain}")
        require("".join(AA[residues[(chain, n)]] for n in numbers) == sequence, f"Wrong sequence {chain}")
        require(all((chain, n, a) in atoms for n in numbers for a in ("N", "CA", "C", "O")), "Missing backbone atom")
        require((chain, len(sequence), "OXT") in atoms, f"Missing free C-terminal OXT {chain}")
    return {key: value for key, value in atoms.items() if key[2] in ("N", "CA", "C", "O", "OXT")}


def displacement(a: dict, b: dict) -> tuple[float, float, int]:
    require(a.keys() == b.keys(), "Backbone atom identity mismatch")
    squared = [sum((a[key][i] - b[key][i]) ** 2 for i in range(3)) for key in a]
    return math.sqrt(statistics.mean(squared)), math.sqrt(max(squared)), len(squared)


def validate_protocol(protocol: dict) -> None:
    require(protocol["state"] == "FROZEN_BEFORE_DIAGNOSTIC_SCORE", "Protocol not frozen before scoring")
    s = protocol["scientific_protocol"]
    fixed = dict(receptor_sequence=RECEPTOR, parents=PARENTS, mutations=MUTATIONS, seeds=list(SEEDS), models=list(MODELS),
                 numberOfRuns=5, moveNeighbours=True, temperature_K=298, pH=7.4, ionStrength_M=.15, water="-IGNORE",
                 noNterm=False, noCterm=False, chains=["A", "B"], expected_jobs=24, expected_pairs=120, expected_analysecomplex_outputs=240,
                 metric="AnalyseComplex A/B Interaction Energy", gain="WT reference energy - mutant energy", units="kcal/mol",
                 aggregation="median five matched gains within each source backbone, then median three complete backbone medians per seed",
                 diagnostic_direction_rule="all three complete backbone medians per mutation/seed must have the experimentally expected sign; mixed or zero is inconclusive",
                 role="POST_SCREEN_DEVELOPMENT_DIAGNOSTIC; not holdout, not experimental free-energy validation")
    require(all(s.get(key) == value for key, value in fixed.items()), "Frozen mutation/chemistry/scoring protocol mismatch")
    for m in s["mutations"]:
        parent = s["parents"][m["parent_id"]]
        position = m["position"]
        require(parent[position - 1] == m["from"] and parent[:position - 1] + m["to"] + parent[position:] == m["sequence"], "Invalid mutation identity")
        require(m["individual_list"] == f"{m['from']}B{position}{m['to']};", "Malformed mutation specification")


def prepare(output: Path = OUTPUT, root: Path = ROOT) -> dict:
    """Freeze inputs and implementation; never edit the already frozen protocol."""
    output, root = output.resolve(), root.resolve()
    protocol_path = output / "protocol.json"
    require(protocol_path.is_file(), "Root must supply frozen protocol.json before preparation")
    protocol = load_json(protocol_path)
    validate_protocol(protocol)
    check_hashes(protocol["source_hashes"], root)
    rows = []
    for mutation in MUTATIONS:
        for seed in SEEDS:
            for model in MODELS:
                job_id = f"{mutation['peptide_id']}_seed{seed}_model{model}"
                source = Path(protocol["source_root"]) / "jobs" / f"{mutation['parent_id']}_seed{seed}_model{model}" / "input_Repair.pdb"
                source_name = source.as_posix()
                require(source_name in protocol["source_hashes"], "Source backbone lacks frozen hash")
                data = (root / source).read_bytes()
                validate_pdb(data.decode("ascii"), PARENTS[mutation["parent_id"]])
                row = dict(job_id=job_id, peptide_id=mutation["peptide_id"], parent_id=mutation["parent_id"], sequence=mutation["sequence"],
                           parent_sequence=PARENTS[mutation["parent_id"]], individual_list=mutation["individual_list"], experimental_direction=mutation["experimental_direction"],
                           seed=seed, model_number=model, source_path=source_name, source_sha256=sha(root / source), input_path=f"jobs/{job_id}/input.pdb",
                           input_sha256=hashlib.sha256(data).hexdigest(), mutation_sha256=hashlib.sha256((mutation["individual_list"] + "\n").encode()).hexdigest(),
                           build_config_sha256=hashlib.sha256(BUILD_CONFIG.encode()).hexdigest(), analyse_template_sha256=hashlib.sha256(ANALYSE_TEMPLATE.encode()).hexdigest())
                rows.append(row)
    manifest = csv_text(rows, MANIFEST_FIELDS)
    preparation_path = output / "preparation.json"
    if preparation_path.exists():
        preparation = load_json(preparation_path)
        require(preparation["protocol_sha256"] == sha(protocol_path), "Frozen protocol changed")
        require(preparation["implementation_hashes"] == implementation_hashes(), "Frozen implementation changed; refusing resume")
        check_hashes(preparation["prepared_hashes"], output)
        require((output / "input_manifest.csv").read_text(encoding="utf-8") == manifest, "Input manifest disagrees with frozen sources")
        return preparation
    require(not (output / "jobs").exists() and not (output / "input_manifest.csv").exists(), "Untracked prepared jobs; refusing overwrite")
    prepared_hashes = {}
    for row in rows:
        job = output / "jobs" / row["job_id"]
        job.mkdir(parents=True)
        shutil.copyfile(root / row["source_path"], job / "input.pdb")
        for name, text in (("individual_list.txt", row["individual_list"] + "\n"), ("build.cfg", BUILD_CONFIG), ("analyse_template.cfg", ANALYSE_TEMPLATE)):
            (job / name).write_text(text, encoding="ascii", newline="\n")
        prepared_hashes.update({(p.relative_to(output)).as_posix(): sha(p) for p in job.iterdir()})
    (output / "input_manifest.csv").write_text(manifest, encoding="utf-8", newline="\n")
    prepared_hashes["input_manifest.csv"] = sha(output / "input_manifest.csv")
    preparation = dict(protocol_sha256=sha(protocol_path), implementation_hashes=implementation_hashes(), prepared_hashes=prepared_hashes,
                       pair_mapping="Observed official PdbList: input_1_0..4.pdb and matching explicitly listed WT_input_1_0..4.pdb", jobs=24)
    save_json(preparation_path, preparation)
    save_json(output / "execution_state.json", dict(state="prepared", preparation_sha256=sha(preparation_path), jobs={}))
    return preparation


def parse_pdblist(path: Path, job: Path, mutation: dict) -> list[dict]:
    require(path.name == "PdbList_input.fxout", "Wrong BuildModel PdbList")
    names = [line.strip() for line in path.read_text(encoding="ascii").splitlines() if line.strip()]
    expected_mut = [f"input_1_{i}.pdb" for i in range(5)]
    expected = set(expected_mut + ["WT_" + name for name in expected_mut])
    require(len(names) == 10 and len(set(names)) == 10 and set(names) == expected, "Missing/ambiguous/unlisted PdbList pair")
    parent_id = mutation["parent_id"]
    original = validate_pdb((job / "input.pdb").read_text(encoding="ascii"), PARENTS[parent_id])
    pairs = []
    for i, name in enumerate(expected_mut, 1):
        mutant_path, wt_path = job / name, job / ("WT_" + name)
        require(mutant_path.is_file() and wt_path.is_file(), "Listed PDB missing")
        wt = validate_pdb(wt_path.read_text(encoding="ascii"), PARENTS[parent_id])
        mutant = validate_pdb(mutant_path.read_text(encoding="ascii"), mutation["sequence"])
        record = dict(repeat=i, wt_path=wt_path.name, mutant_path=mutant_path.name)
        for prefix, a, b in (("wt_source", original, wt), ("mutant_source", original, mutant), ("pair", wt, mutant)):
            rmsd, maximum, n = displacement(a, b)
            record.update({f"{prefix}_backbone_rmsd_A": rmsd, f"{prefix}_backbone_max_A": maximum, "backbone_atom_count": n})
        pairs.append(record)
    return pairs


def parse_interaction_terms(path: Path, expected_pdb: str) -> dict[str, float]:
    require(path.name.startswith(("Interaction_", "Summary_")) and path.suffix == ".fxout", "Not an AnalyseComplex interaction table")
    header, result = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        cells = [c.strip() for c in line.split("\t")]
        if "Interaction Energy" in cells:
            require(header is None and len(cells) == len(set(cells)), "Ambiguous duplicate energy header")
            require(all(k in cells for k in ("Pdb", "Group1", "Group2")), "Missing PDB/chain headers")
            header = cells
            continue
        if header is None or not line.strip():
            continue
        require(len(cells) == len(header), "Malformed energy row")
        row = dict(zip(header, cells))
        require({row["Group1"], row["Group2"]} == {"A", "B"}, "Wrong interaction chains")
        require(row["Pdb"] in (expected_pdb, "./" + expected_pdb), "Wrong scored PDB")
        terms = {k: float(v) for k, v in row.items() if k not in ("Pdb", "Group1", "Group2")}
        require(all(math.isfinite(v) for v in terms.values()), "Nonfinite energy term")
        result.append(terms)
    require(header is not None and len(result) == 1, "Expected exactly one AnalyseComplex A/B row")
    return result[0]


def raw_energy(directory: Path, name: str) -> dict:
    interaction, summary = sorted(directory.glob("Interaction_*.fxout")), sorted(directory.glob("Summary_*.fxout"))
    require(len(interaction) == 1 and len(summary) == 1, "Missing/ambiguous actual AnalyseComplex Interaction/Summary outputs")
    terms = parse_interaction_terms(interaction[0], name)
    require("IntraclashesGroup1" in terms and "IntraclashesGroup2" in terms, "Missing intraclashes terms")
    tables = {interaction[0].name: terms}
    for path in summary:
        extra = parse_interaction_terms(path, name)
        require(math.isclose(extra["Interaction Energy"], terms["Interaction Energy"], rel_tol=1e-8, abs_tol=1e-6), "Summary/interaction energy disagreement")
        tables[path.name] = extra
    return dict(energy_kcal_mol=terms["Interaction Energy"], tables=tables)


def matched_gain(wt: float, mutant: float) -> float:
    require(math.isfinite(wt) and math.isfinite(mutant), "Nonfinite paired energy")
    return wt - mutant


def summarize_backbone(gains: list[float]) -> dict:
    require(all(math.isfinite(v) for v in gains), "Nonfinite gain")
    require(len(gains) <= 5, "Too many repeat pairs")
    complete = len(gains) == 5
    return dict(status="complete" if complete else "incomplete", complete=complete, n_complete_pairs=len(gains),
                median_gain_kcal_mol=statistics.median(gains) if complete else "", minimum_gain_kcal_mol=min(gains) if gains else "",
                maximum_gain_kcal_mol=max(gains) if gains else "", gain_range_kcal_mol=max(gains) - min(gains) if gains else "",
                positive_pairs=sum(v > 0 for v in gains), negative_pairs=sum(v < 0 for v in gains), zero_pairs=sum(v == 0 for v in gains))


def summarize_seed(medians: list[float], expected_direction: str) -> dict:
    require(expected_direction in ("stronger", "weaker") and len(medians) <= 3 and all(math.isfinite(v) for v in medians), "Invalid seed summary")
    expected = 1 if expected_direction == "stronger" else -1
    complete = len(medians) == 3
    label = "not_evaluable"
    if complete:
        label = "correct" if all(v * expected > 0 for v in medians) else "incorrect" if all(v * expected < 0 for v in medians) else "inconclusive"
    return dict(complete=complete, n_complete_backbones=len(medians), median_gain_kcal_mol=statistics.median(medians) if complete else "",
                minimum_backbone_gain_kcal_mol=min(medians) if medians else "", maximum_backbone_gain_kcal_mol=max(medians) if medians else "",
                positive_backbones=sum(v > 0 for v in medians), negative_backbones=sum(v < 0 for v in medians), zero_backbones=sum(v == 0 for v in medians), diagnostic_direction=label)


def command_log(command: list[str], cwd: Path, stage: str, timeout: float) -> dict:
    start = time.monotonic()
    metadata = dict(command=command, cwd=str(cwd), started_utc=datetime.now(timezone.utc).isoformat())
    try:
        process = subprocess.run(command, cwd=cwd, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout,
                                 env={**os.environ, "LC_ALL": "C", "LANG": "C"}, shell=False)
        stdout, stderr = process.stdout, process.stderr
        metadata.update(returncode=process.returncode, error="")
    except (OSError, subprocess.TimeoutExpired) as error:
        stdout, stderr = getattr(error, "stdout", "") or "", getattr(error, "stderr", "") or ""
        stdout = stdout.decode(errors="replace") if isinstance(stdout, bytes) else stdout
        stderr = stderr.decode(errors="replace") if isinstance(stderr, bytes) else stderr
        metadata.update(returncode=None, error=str(error))
    (cwd / f"{stage}.stdout.txt").write_text(stdout, encoding="utf-8")
    (cwd / f"{stage}.stderr.txt").write_text(stderr, encoding="utf-8")
    metadata["seconds"] = time.monotonic() - start
    save_json(cwd / f"{stage}.command.json", metadata)
    return metadata


def job_artifacts(job: Path) -> dict:
    return {p.relative_to(job).as_posix(): sha(p) for p in sorted(job.rglob("*")) if p.is_file() and p.relative_to(job).as_posix() not in PREPARED_NAMES | {"result.json"}}


def archive_interrupted(job: Path) -> None:
    leftovers = [p for p in job.iterdir() if p.name not in PREPARED_NAMES | {"attempts", "result.json"}]
    if not leftovers:
        return
    archive = job / "attempts" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
    archive.mkdir(parents=True)
    for path in leftovers:
        shutil.move(str(path), str(archive / path.name))
    save_json(archive / "artifact_hashes.json", {p.relative_to(archive).as_posix(): sha(p) for p in sorted(archive.rglob("*")) if p.is_file()})


def role_directory(job: Path, repeat: int, role: str) -> Path:
    return job / "ac" / f"repeat{repeat}" / role


def validate_job(job: Path, row: dict, signature: str, expected_hash: str) -> dict:
    require(sha(job / "result.json") == expected_hash, "Tracked job result changed")
    result = load_json(job / "result.json")
    require(result["signature"] == signature and result["job_id"] == row["job_id"], "Wrong result identity/signature")
    require(result["artifact_hashes"] == job_artifacts(job), "Raw/artifact hashes changed")
    require(result["status"] in ("ok", "failed") and len(result["pairs"]) == 5, "Invalid job completeness")
    require(type(result["provenance_valid"]) is bool, "Invalid job provenance state")
    if not result["provenance_valid"]:
        require(result["status"] == "failed" and result["failure_reason"] and all(p["status"] == "failed" and "gain_kcal_mol" not in p for p in result["pairs"]),
                "Provenance failure cannot carry successful pairs or gains")
    require([p["repeat"] for p in result["pairs"]] == [1, 2, 3, 4, 5], "Ambiguous cached repeat identity")
    if result["build_valid"]:
        fresh = parse_pdblist(job / "PdbList_input.fxout", job, row)
        require(load_json(job / "build.command.json")["returncode"] == 0, "BuildModel did not succeed")
        for pair, mapping in zip(result["pairs"], fresh):
            require(all(pair.get(k) == v for k, v in mapping.items()), "Cached PDB mapping/displacement differs from raw")
            require(set(pair["roles"]) == {"wt_reference", "mutant"}, "Missing/ambiguous pair role identities")
            for role, name in (("wt_reference", mapping["wt_path"]), ("mutant", mapping["mutant_path"])):
                role_result = pair["roles"][role]
                if role_result["status"] == "ok":
                    directory = role_directory(job, pair["repeat"], role)
                    require(sha(directory / name) == sha(job / name), "Scored PDB differs from specific BuildModel pair")
                    require(load_json(directory / "analyse.command.json")["returncode"] == 0, "AnalyseComplex did not succeed")
                    require(role_result["raw"] == raw_energy(directory, name), "Cached energy terms differ from raw")
                else:
                    require(role_result["status"] in ("failed", "not_run") and "raw" not in role_result, "Failure carries energy")
            both = all(v["status"] == "ok" for v in pair["roles"].values())
            require((pair["status"] == "ok") == (both and result["provenance_valid"]), "Pair success disagrees with role results/provenance")
            if pair["status"] == "ok":
                require(pair["gain_kcal_mol"] == matched_gain(pair["roles"]["wt_reference"]["raw"]["energy_kcal_mol"], pair["roles"]["mutant"]["raw"]["energy_kcal_mol"]), "Cached paired gain differs from raw")
            else:
                require(pair["status"] == "failed" and pair["failure_reason"] and "gain_kcal_mol" not in pair, "Failed pair must have an explicit reason and no gain")
    else:
        require(all(p["status"] == "failed" and "gain_kcal_mol" not in p for p in result["pairs"]), "Invalid failed BuildModel result")
    require((result["status"] == "ok") == all(p["status"] == "ok" for p in result["pairs"]), "Job marked complete with incomplete pairs")
    return result


def run(output: Path = OUTPUT, foldx: Path | None = None, root: Path = ROOT, timeout: float = 1800, workers: int = 2, runner=command_log) -> dict:
    """At most two independent jobs; each BuildModel and every AC has isolated cwd."""
    preparation = prepare(output, root)
    output, root = output.resolve(), root.resolve()
    require(workers in (1, 2), "Only one or two isolated CPU jobs permitted")
    require(foldx is not None and foldx.is_file() and foldx.is_absolute(), "Existing absolute licensed FoldX executable required")
    binary = foldx.resolve()
    protocol = load_json(output / "protocol.json")
    try:
        binary_relative = binary.relative_to(root).as_posix()
    except ValueError as error:
        raise ValueError("FoldX executable is outside frozen source root") from error
    require(protocol["source_hashes"].get(binary_relative) == sha(binary), "Executable not frozen by protocol")
    resources = {str(Path(name).relative_to(binary.parent.relative_to(root))).replace("\\", "/"): name for name in protocol["source_hashes"]
                 if (root / name) != binary and (root / name).is_relative_to(binary.parent)}
    execution = dict(binary_relative=binary_relative, binary_sha256=sha(binary), resources=resources,
                     preparation_sha256=sha(output / "preparation.json"), workers_max=2)
    execution_path = output / "execution_manifest.json"
    if execution_path.exists():
        require(load_json(execution_path) == execution, "Frozen execution settings changed")
    else:
        save_json(execution_path, execution)
    signature = hashlib.sha256(json_text(execution).encode()).hexdigest()
    state = load_json(output / "execution_state.json")
    require(state["preparation_sha256"] == sha(output / "preparation.json"), "Preparation manifest changed")
    rows = read_csv(output / "input_manifest.csv")
    pending = []
    for row in rows:
        job = output / "jobs" / row["job_id"]
        if row["job_id"] in state["jobs"]:
            validate_job(job, row, signature, state["jobs"][row["job_id"]])
        else:
            require(not (job / "result.json").exists(), "Untracked result; refusing overwrite")
            pending.append(row)
    state["state"] = "running"
    save_json(output / "execution_state.json", state)
    lock = threading.Lock()

    def copy_resources(directory: Path) -> None:
        for relative, source in resources.items():
            target = directory / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(root / source, target)

    def score_job(row: dict) -> None:
        job = output / "jobs" / row["job_id"]
        archive_interrupted(job)
        result = dict(job_id=row["job_id"], signature=signature, status="failed", build_valid=False, provenance_valid=True, failure_reason="",
                      pairs=[dict(repeat=i, status="failed", failure_reason="Job did not finish scoring", roles={role: dict(status="not_run", failure_reason="Not attempted") for role in ("wt_reference", "mutant")}) for i in range(1, 6)])
        phase = "preflight_provenance"
        try:
            check_hashes(protocol["source_hashes"], root)
            check_hashes(preparation["prepared_hashes"], output)
            phase = "build"
            copy_resources(job)
            metadata = runner([str(binary), "-f", "build.cfg"], job, "build", timeout)
            require(metadata["returncode"] == 0, f"BuildModel failed: {metadata.get('error') or metadata['returncode']}")
            pairs = parse_pdblist(job / "PdbList_input.fxout", job, row)
            result["build_valid"] = True
            for cached, mapping in zip(result["pairs"], pairs):
                cached.update(mapping)
            phase = "scoring"
            for pair in result["pairs"]:
                for role, name in (("wt_reference", pair["wt_path"]), ("mutant", pair["mutant_path"])):
                    directory = role_directory(job, pair["repeat"], role)
                    role_result = dict(status="failed", failure_reason="")
                    try:
                        directory.mkdir(parents=True)
                        shutil.copyfile(job / name, directory / name)
                        (directory / "analyse.cfg").write_text(ANALYSE_TEMPLATE.format(pdb=name), encoding="ascii", newline="\n")
                        copy_resources(directory)
                        metadata = runner([str(binary), "-f", "analyse.cfg"], directory, "analyse", timeout)
                        require(metadata["returncode"] == 0, f"AnalyseComplex failed: {metadata.get('error') or metadata['returncode']}")
                        require(sha(directory / name) == sha(job / name), "Scored PDB changed during AnalyseComplex")
                        role_result.update(status="ok", raw=raw_energy(directory, name))
                    except (ValueError, OSError) as error:
                        role_result["failure_reason"] = str(error)
                    pair["roles"][role] = role_result
                if all(v["status"] == "ok" for v in pair["roles"].values()):
                    pair.update(status="ok", failure_reason="", gain_kcal_mol=matched_gain(pair["roles"]["wt_reference"]["raw"]["energy_kcal_mol"], pair["roles"]["mutant"]["raw"]["energy_kcal_mol"]))
                else:
                    pair["failure_reason"] = "; ".join(f"{role}: {v['failure_reason']}" for role, v in pair["roles"].items() if v["status"] != "ok")
            phase = "postflight_provenance"
            check_hashes(protocol["source_hashes"], root)
            check_hashes(preparation["prepared_hashes"], output)
            for directory in [job] + [role_directory(job, p["repeat"], role) for p in result["pairs"] for role in ("wt_reference", "mutant")]:
                for relative, source in resources.items():
                    require(sha(directory / relative) == protocol["source_hashes"][source], "Frozen tool resource changed during command")
            result["status"] = "ok" if all(p["status"] == "ok" for p in result["pairs"]) else "failed"
            if result["status"] == "failed":
                result["failure_reason"] = "; ".join(f"repeat {p['repeat']}: {p['failure_reason']}" for p in result["pairs"] if p["status"] != "ok")
        except (ValueError, OSError) as error:
            result["failure_reason"] = str(error)
            if phase in ("preflight_provenance", "postflight_provenance"):
                result["provenance_valid"] = False
            for pair in result["pairs"]:
                if not result["provenance_valid"] or pair["status"] != "ok":
                    pair.update(status="failed", failure_reason=f"{phase}: {error}")
                    pair.pop("gain_kcal_mol", None)
                    for role_result in pair["roles"].values():
                        if role_result["status"] == "not_run":
                            role_result["failure_reason"] = f"Job stopped during {phase}: {error}"
        result["artifact_hashes"] = job_artifacts(job)
        save_json(job / "result.json", result)
        with lock:
            state["jobs"][row["job_id"]] = sha(job / "result.json")
            save_json(output / "execution_state.json", state)
        print(f"{row['job_id']}: {result['status']}", flush=True)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        list(pool.map(score_job, pending))
    results = [load_json(output / "jobs" / row["job_id"] / "result.json") for row in rows]
    state.update(state="complete" if all(r["status"] == "ok" for r in results) else "complete_with_failures",
                 job_counts={status: sum(r["status"] == status for r in results) for status in ("ok", "failed")})
    save_json(output / "execution_state.json", state)
    return state


def old_energies(protocol: dict, root: Path) -> dict:
    rows = read_csv(root / protocol["prior_result_csv"])
    result = {}
    for row in rows:
        key = (row["peptide_id"], int(row["seed"]), int(row["model_number"]))
        require(key not in result, "Ambiguous old independent AFM model identity")
        value = float(row["interaction_energy_kcal_mol"]) if row["status"] == "ok" else None
        require(value is None or math.isfinite(value), "Nonfinite old independent energy")
        result[key] = dict(energy=value, quality_pass=row.get("quality_pass", ""), sequence=row["sequence"])
    for mutation in MUTATIONS:
        for pid, seq in ((mutation["parent_id"], PARENTS[mutation["parent_id"]]), (mutation["peptide_id"], mutation["sequence"])):
            for seed in SEEDS:
                for model in MODELS:
                    require((pid, seed, model) in result and result[(pid, seed, model)]["sequence"] == seq, "Missing/wrong old comparison identity")
    return result


def analysis_tables(output: Path, root: Path = ROOT) -> tuple[dict, dict]:
    preparation = prepare(output, root)
    protocol = load_json(output / "protocol.json")
    state = load_json(output / "execution_state.json")
    require(state["preparation_sha256"] == sha(output / "preparation.json"), "Frozen preparation changed")
    execution = load_json(output / "execution_manifest.json") if (output / "execution_manifest.json").exists() else None
    signature = hashlib.sha256(json_text(execution).encode()).hexdigest() if execution else None
    if execution:
        require(execution["preparation_sha256"] == sha(output / "preparation.json"), "Execution manifest preparation mismatch")
        require(protocol["source_hashes"].get(execution["binary_relative"]) == execution["binary_sha256"], "Execution executable fingerprint mismatch")
    old = old_energies(protocol, root)
    pairs, terms, backbones, jobs = [], [], [], []
    ac_counts = {"ok": 0, "failed": 0, "not_run": 0}
    manifest = read_csv(output / "input_manifest.csv")
    require(set(state["jobs"]).issubset({r["job_id"] for r in manifest}), "Unrecognized tracked jobs")
    for row in manifest:
        job = output / "jobs" / row["job_id"]
        result = None
        if row["job_id"] in state["jobs"]:
            require(execution is not None, "Results lack frozen execution manifest")
            result = validate_job(job, row, signature, state["jobs"][row["job_id"]])
        else:
            require(not (job / "result.json").exists(), "Untracked job result")
        job_status = result["status"] if result else "not_run"
        jobs.append(job_status)
        old_parent = old[(row["parent_id"], int(row["seed"]), int(row["model_number"]))]["energy"]
        old_mutant = old[(row["peptide_id"], int(row["seed"]), int(row["model_number"]))]["energy"]
        old_gain = matched_gain(old_parent, old_mutant) if old_parent is not None and old_mutant is not None else ""
        group = []
        for repeat in range(1, 6):
            cached = result["pairs"][repeat - 1] if result else dict(repeat=repeat, status="not_run", failure_reason="No real FoldX result", roles={})
            record = {key: "" for key in PAIR_FIELDS}
            record.update({key: row[key] for key in ("job_id", "peptide_id", "parent_id", "seed", "model_number", "experimental_direction", "source_path")})
            record.update(repeat=repeat, status=cached["status"], failure_reason=cached["failure_reason"], old_parent_energy_kcal_mol=old_parent if old_parent is not None else "",
                          old_mutant_energy_kcal_mol=old_mutant if old_mutant is not None else "", old_independent_gain_kcal_mol=old_gain,
                          provenance_valid=result["provenance_valid"] if result else "", job_failure_reason=result["failure_reason"] if result else "", units="kcal/mol")
            record.update({key: value for key, value in cached.items() if key in PAIR_FIELDS})
            if cached.get("wt_path"):
                record["wt_path"] = (job / cached["wt_path"]).relative_to(output).as_posix()
                record["mutant_path"] = (job / cached["mutant_path"]).relative_to(output).as_posix()
            for role, prefix in (("wt_reference", "wt"), ("mutant", "mutant")):
                r = cached.get("roles", {}).get(role, dict(status="not_run"))
                ac_counts[r["status"]] += 1
                if r["status"] == "ok":
                    record[f"{prefix}_energy_kcal_mol"] = r["raw"]["energy_kcal_mol"]
                    raw_paths = []
                    for filename, table_terms in r["raw"]["tables"].items():
                        path = (role_directory(job, repeat, role) / filename).relative_to(output).as_posix()
                        raw_paths.append(path)
                        for term, value in table_terms.items():
                            terms.append(dict(job_id=row["job_id"], peptide_id=row["peptide_id"], parent_id=row["parent_id"], seed=row["seed"], model_number=row["model_number"],
                                              repeat=repeat, role=role, table_path=path, term=term, value=value, units="count" if "Residues" in term else "kcal/mol"))
                    record[f"{prefix}_raw_paths"] = ";".join(raw_paths)
                    comparison = old_parent if role == "wt_reference" else old_mutant
                    if comparison is not None and result["provenance_valid"]:
                        record["wt_repacking_delta_kcal_mol" if role == "wt_reference" else "mutant_vs_old_mutant_delta_kcal_mol"] = r["raw"]["energy_kcal_mol"] - comparison
            pairs.append(record)
            group.append(record)
        gains = [p["gain_kcal_mol"] for p in group if p["status"] == "ok"]
        summary = {key: row[key] for key in ("job_id", "peptide_id", "parent_id", "seed", "model_number", "experimental_direction")}
        summary.update(summarize_backbone(gains), job_status=job_status, failed_pairs=sum(p["status"] == "failed" for p in group),
                       provenance_valid=result["provenance_valid"] if result else "", job_failure_reason=result["failure_reason"] if result else "",
                       not_run_pairs=sum(p["status"] == "not_run" for p in group), old_independent_gain_kcal_mol=old_gain,
                       old_parent_quality_pass=old[(row["parent_id"], int(row["seed"]), int(row["model_number"]))]["quality_pass"],
                       old_mutant_quality_pass=old[(row["peptide_id"], int(row["seed"]), int(row["model_number"]))]["quality_pass"], units="kcal/mol")
        for prefix in ("wt", "mutant"):
            values = [p[f"{prefix}_energy_kcal_mol"] for p in group if p[f"{prefix}_energy_kcal_mol"] != ""]
            summary[f"{prefix}_energy_range_kcal_mol"] = max(values) - min(values) if values else ""
            summary[f"{prefix}_energy_median_kcal_mol"] = statistics.median(values) if len(values) == 5 else ""
        deltas = [p["wt_repacking_delta_kcal_mol"] for p in group if p["wt_repacking_delta_kcal_mol"] != ""]
        summary["wt_repacking_delta_median_kcal_mol"] = statistics.median(deltas) if len(deltas) == 5 else ""
        summary["matched_minus_old_gain_kcal_mol"] = summary["median_gain_kcal_mol"] - old_gain if summary["complete"] and old_gain != "" else ""
        backbones.append(summary)
    seeds = []
    for mutation in MUTATIONS:
        for seed in SEEDS:
            selected = [b for b in backbones if b["peptide_id"] == mutation["peptide_id"] and int(b["seed"]) == seed]
            medians = [b["median_gain_kcal_mol"] for b in selected if b["complete"]]
            record = dict(peptide_id=mutation["peptide_id"], parent_id=mutation["parent_id"], seed=seed, experimental_direction=mutation["experimental_direction"])
            record.update(summarize_seed(medians, mutation["experimental_direction"]), complete_models=";".join(str(b["model_number"]) for b in selected if b["complete"]),
                          failed_pairs=sum(b["failed_pairs"] for b in selected), not_run_pairs=sum(b["not_run_pairs"] for b in selected), units="kcal/mol")
            old_gains = [b["old_independent_gain_kcal_mol"] for b in selected if b["old_independent_gain_kcal_mol"] != ""]
            old_summary = summarize_seed(old_gains, mutation["experimental_direction"])
            record.update(old_independent_median_gain_kcal_mol=old_summary["median_gain_kcal_mol"], old_independent_direction=old_summary["diagnostic_direction"],
                          matched_minus_old_median_gain_kcal_mol=record["median_gain_kcal_mol"] - old_summary["median_gain_kcal_mol"] if record["complete"] and old_summary["complete"] else "")
            seeds.append(record)
    counts = {status: sum(p["status"] == status for p in pairs) for status in ("ok", "failed", "not_run")}
    job_counts = {status: jobs.count(status) for status in ("ok", "failed", "not_run")}
    complete = job_counts["ok"] == 24 and counts["ok"] == 120 and ac_counts["ok"] == 240
    decision = dict(role=protocol["scientific_protocol"]["role"], complete=complete, job_counts=job_counts, pair_counts=counts,
                    provenance_failed_jobs=sum(b["provenance_valid"] is False for b in backbones),
                    analysecomplex_counts=ac_counts, expected_jobs=24, expected_pairs=120, expected_analysecomplex_outputs=240,
                    real_scoring="COMPLETE" if complete else "NOT_RUN" if job_counts["not_run"] == 24 else "PARTIAL" if job_counts["not_run"] else "COMPLETE_WITH_FAILURES",
                    diagnostic_direction_counts={label: sum(s["diagnostic_direction"] == label for s in seeds) for label in ("correct", "incorrect", "inconclusive", "not_evaluable")},
                    seed_diagnostics=seeds, replaces_old_screen_gate=False, experimental_Kd="NOT_MEASURED", holdout_validation="NOT_RUN",
                    microscopic_gain_caution="Signs close to printed precision or FoldX repeat variation are descriptive, even when all three backbone signs agree; they do not establish validated affinity improvement",
                    backbone_displacement="Measured directly in original coordinate frame for N/CA/C/O/OXT; no alignment or assumption of fixed backbone",
                    gain="WT reference Interaction Energy - mutant Interaction Energy", units="kcal/mol",
                    limitations=protocol["scientific_protocol"]["limitations"])
    return {"pair_results.csv": pairs, "backbone_summary.csv": backbones, "seed_summary.csv": seeds, "energy_terms.csv": terms}, decision


def chinese_report(decision: dict) -> str:
    lines = ["# FoldX matched WT/mutant DEVELOPMENT 诊断", "", "本诊断来自已看过结果的开发面板；不是 holdout 或新实验 Kd 验证，不能替换旧筛选门槛。", "",
             f"真实评分状态：{decision['real_scoring']}；任务={decision['job_counts']}；配对={decision['pair_counts']}；AnalyseComplex={decision['analysecomplex_counts']}；来源/资源完整性失败任务={decision['provenance_failed_jobs']}；完整={decision['complete']}。", "",
             "来源/资源完整性校验失败时，保留原始输出用于审计，五个配对全部标记失败并禁止产生 gain；随后恢复冻结文件也不会把已失败记录提升为成功。", "",
             "每个 mutant 仅与实际 PdbList 明确列出的同后缀 WT reference 配对。gain = E(WT reference) − E(mutant)，正值表示更有利；能量仅来自 AnalyseComplex A/B Interaction Energy（kcal/mol）。BuildModel Dif Total Energy 不作为结合能。", "",
             "每个原始 backbone 要求五个完整 gain 后取中位数；每个 seed 要求三个完整 backbone 中位数后取中位数。三个 backbone 均符合实验方向才记 correct；全相反为 incorrect；混合或零为 inconclusive；不完整为 not_evaluable。失败/缺失不选择性剔除。", ""]
    for seed in decision["seed_diagnostics"]:
        lines.append(f"{seed['peptide_id']} seed {seed['seed']}：{seed['diagnostic_direction']}；完整 backbone {seed['n_complete_backbones']}/3；matched gain={seed['median_gain_kcal_mol']}；旧独立 AFM gain={seed['old_independent_median_gain_kcal_mol']} ({seed['old_independent_direction']})。")
    lines += ["", "CSV 保留各重复配对、能量范围、原始数值项（含两个 intraclashes）、旧 AFM 同 mutation/parent/seed/model 差值，以及 WT reference 相对旧 WT 的 repacking 能量变化。骨架 N/CA/C/O/OXT 相对来源和配对内的位移均直接测量；不预设为零。", "",
              "共同起始结构与局部优化限制了构象变化，不能充分模拟螺旋重排或游离肽构象与结合前后的全自由能，尤其可能漏掉 P27A/N8A 的构象机制。五个 FoldX 重复和三个 AFM 参数模型均不是生物学重复。方向改善只能提示几何/局部 repacking 敏感性，不能证明 AFM 因果机制或可靠亲和力预测。", "",
              "接近输出数值舍入精度或重复波动的微小 gain，即使三个 backbone 同号，也仅是描述性方向；不能据此声称可靠亲和力改善。保留实际打印数值与严格零值规则，不另设事后阈值。", ""]
    return "\n".join(lines)


def published_contents(tables: dict, decision: dict) -> dict[str, str]:
    contents = {name: csv_text(rows, PAIR_FIELDS if name == "pair_results.csv" else TERM_FIELDS if name == "energy_terms.csv" else None) for name, rows in tables.items()}
    contents.update({"decision.json": json_text(decision), "report_zh.md": chinese_report(decision)})
    return contents


def analyze(output: Path = OUTPUT, root: Path = ROOT) -> dict:
    tables, decision = analysis_tables(output, root)
    directory = output / "analysis"
    directory.mkdir(exist_ok=True)
    for name, text in published_contents(tables, decision).items():
        (directory / name).write_text(text, encoding="utf-8", newline="\n")
    save_json(directory / "artifact_hashes.json", {name: sha(directory / name) for name in published_contents(tables, decision)})
    return decision


def verify(output: Path = OUTPUT, root: Path = ROOT) -> dict:
    tables, decision = analysis_tables(output, root)
    directory = output / "analysis"
    published = (directory / "decision.json").exists()
    if published:
        check_hashes(load_json(directory / "artifact_hashes.json"), directory)
        for name, text in published_contents(tables, decision).items():
            require((directory / name).read_text(encoding="utf-8") == text, f"Published artifact differs from raw-recomputed values: {name}")
    result = dict(inputs="PASS", frozen_hashes="PASS", raw_result_reparse="PASS" if decision["analysecomplex_counts"]["ok"] else "NOT_APPLICABLE_NO_SCORED_PAIRS",
                  published_tables="PASS" if published else "NOT_PUBLISHED", **decision)
    save_json(output / "verification.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "run", "analyze", "verify"))
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--foldx", type=Path)
    parser.add_argument("--timeout", type=float, default=1800)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    try:
        if args.command == "run":
            result = run(args.output.resolve(), args.foldx.resolve() if args.foldx else None, args.root.resolve(), args.timeout, args.workers)
        else:
            result = globals()[args.command](args.output.resolve(), args.root.resolve())
        print(json_text(result))
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"Matched FoldX diagnostic stopped: {error}\n")


if __name__ == "__main__":
    main()
