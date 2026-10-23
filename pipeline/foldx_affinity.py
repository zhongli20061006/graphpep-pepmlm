"""Frozen DEVELOPMENT/SCREEN MDM2 panel; only real AnalyseComplex output is scored."""
from __future__ import annotations

import argparse
import csv
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
import time
from datetime import datetime, timezone

from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "output/affinity_diagnostic_mdm2_20261006"
OUTPUT = ROOT / "output/foldx_affinity_mdm2_20261006"
RECEPTOR = "ETLVRPKPLLLKLLKSVGAQKDTYTMKEVLFYLGQYIMTKRLYDEKQQHIVYCSNDLLGDLFGVPSFSVKEHRKIYTMIYRNLVV"
SEEDS, MODELS = (20261011, 20261012), (1, 2, 3)
REFERENCES = {
    "PMI_WT": ("TSFAEYWNLLSP", 3.2, 1.1, "PMI_WT", 0),
    "PMI_N8A": ("TSFAEYWALLSP", .49, .21, "PMI_WT", 8),
    "PMI_Y6A": ("TSFAEAWNLLSP", 610., 70., "PMI_WT", 6),
    "PMI_F3A": ("TSAAEYWNLLSP", 38000., 2000., "PMI_WT", 3),
    "PMI_W7A": ("TSFAEYANLLSP", 160000., 30000., "PMI_WT", 7),
    "P53_WT": ("ETFSDLWKLLPE", 440., 40., "P53_WT", 0),
    "P53_P27A": ("ETFSDLWKLLAE", 51., 3., "P53_WT", 11),
    "P53_K24A": ("ETFSDLWALLPE", 230., 20., "P53_WT", 8),
    "P53_L22A": ("ETFSDAWKLLPE", 5000., 400., "P53_WT", 6),
    "P53_L26A": ("ETFSDLWKLAPE", 27000., 1000., "P53_WT", 10),
}
STRONG = ("PMI_N8A", "PMI_Y6A", "PMI_F3A", "PMI_W7A", "P53_P27A", "P53_L22A", "P53_L26A")
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
SETTINGS = "temperature=298\npH=7.4\nionStrength=0.15\nwater=-IGNORE\nscreen=false\nnoNterm=false\nnoCterm=false\n"
MANIFEST_FIELDS = "job_id peptide_id family sequence parent_id kd_nM log10_kd_nM seed model_number quality_pass quality_fail_reasons source_path source_sha256 input_path input_sha256 repair_config_sha256 analyse_config_sha256 heavy_atoms removed_hydrogens".split()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def digest(content: str) -> str:
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


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


def heavy_pdb(text: str, expected: dict[str, str], full_termini: bool = True) -> tuple[str, dict]:
    """Keep atom lines byte-for-byte (except newline), rejecting ambiguous chemistry."""
    residues, atoms, kept, hydrogens = {}, set(), [], 0
    for line in text.splitlines():
        require(not line.startswith(("MODEL ", "ENDMDL")), "Multiple/model-delimited PDB unsupported")
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        require(len(line) >= 54, "Truncated atom record")
        require(line.startswith("ATOM  "), "HETATM/nonstandard chemistry rejected")
        require(line[16] == " " and line[26] == " ", "Altloc/insertion code rejected")
        name, resname, chain = line[12:16].strip(), line[17:20].strip(), line[21]
        require(chain in expected and resname in AA, "Unexpected chain or nonstandard residue")
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        element = element or re.sub(r"^[0-9]+", "", name).upper()[:1]
        if element in ("H", "D"):
            hydrogens += 1
            continue
        require(element in ("C", "N", "O", "S"), f"Unsupported heavy element {element}")
        number = int(line[22:26])
        key = (chain, number)
        require(key not in residues or residues[key] == resname, "Conflicting residue identity")
        residues[key] = resname
        atom = (*key, name)
        require(atom not in atoms, "Duplicate heavy atom")
        atoms.add(atom)
        require(all(math.isfinite(float(line[start:start + 8])) for start in (30, 38, 46)), "Nonfinite coordinate")
        line.encode("ascii")
        kept.append(line)
    for chain, sequence in expected.items():
        numbers = sorted(number for c, number in residues if c == chain)
        require(numbers == list(range(1, len(sequence) + 1)), f"Incomplete/noncontiguous chain {chain}")
        require("".join(AA[residues[(chain, n)]] for n in numbers) == sequence, f"Wrong sequence {chain}")
        for number in numbers:
            require(all((chain, number, name) in atoms for name in ("N", "CA", "C", "O")), "Missing backbone atom")
        if full_termini:
            require((chain, len(sequence), "OXT") in atoms, f"Missing free C-terminal OXT in {chain}")
    output = []
    previous = None
    for line in kept:
        if previous is not None and line[21] != previous:
            output.append("TER")
        output.append(line)
        previous = line[21]
    return "\n".join(output + ["TER", "END", ""]), {"heavy_atoms": len(atoms), "removed_hydrogens": hydrogens}


def validated_sources(source: Path, wrapper: Path) -> tuple[list[dict], dict]:
    refs = read_csv(source / "reference_affinity.csv")
    old = load_json(source / "protocol.json")
    require(len(refs) == 10 and {r["peptide_id"] for r in refs} == set(REFERENCES), "Wrong reference panel")
    require(old["receptor_sequence"] == RECEPTOR and old["seeds"] == list(SEEDS) and old["models"] == list(MODELS), "Source protocol identities changed")
    gate = {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10., "minimum_distance_min": 1.5}
    require(old["quality_gate"] == gate and old["minimum_passing_models_per_seed"] == 2 and old["strong_direction_cases"] == list(STRONG), "Source criteria changed")
    ref_map = {r["peptide_id"]: r for r in refs}
    for pid, reference in ref_map.items():
        sequence, kd, sd, parent, position = REFERENCES[pid]
        require(reference["sequence"] == reference["peptide_sequence"] == sequence and reference["parent_id"] == parent, f"Reference sequence/parent mismatch {pid}")
        require(int(reference["position"]) == position and reference["family"] == ("PMI" if pid.startswith("PMI") else "p53"), "Reference mutation/family mismatch")
        require(reference["mutation"] == pid.split("_", 1)[1] and int(reference["table"]) == (1 if pid.startswith("PMI") else 3), "Reference mutation/table mismatch")
        for field, value in (("kd_nM", kd), ("sd_nM", sd), ("kd_M", kd * 1e-9), ("sd_M", sd * 1e-9), ("log10_kd_nM", math.log10(kd))):
            require(math.isclose(float(reference[field]), value, rel_tol=1e-12, abs_tol=1e-15), f"Reference units/value mismatch {pid} {field}")
        require(reference["n_terminal"] == "free amine" and reference["c_terminal"] == "free acid" and reference["receptor_construct"] == "MDM2 25-109", "Reference chemistry mismatch")
        require(reference["source_doi"] == "10.1016/j.jmb.2010.03.005" and reference["source_url"] == "https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/" and reference["assay"] == "SPR competition", "Reference provenance mismatch")
    rows = read_csv(source / "analysis/model_results.csv")
    identities = set()
    for row in rows:
        pid, seed, model = row["peptide_id"], int(row["seed"]), int(row["model_number"])
        identity = (pid, seed, model)
        require(pid in ref_map and seed in SEEDS and model in MODELS and identity not in identities, "Invalid/duplicate model identity")
        identities.add(identity)
        require(row["model_id"] == f"seed{seed}_model{model}" and row["receptor_sequence"] == RECEPTOR and row["scope"] == "full" and row["source"] == "AFM", "Wrong source model identity")
        for field, value in ref_map[pid].items():
            same = row[field] == value
            if not same:
                try:
                    same = math.isclose(float(row[field]), float(value), rel_tol=1e-12, abs_tol=1e-15)
                except ValueError:
                    pass
            require(same, f"Model reference mismatch {pid} {field}")
        path = Path(row["complex_path"])
        require(path.is_absolute() and sha(path) == row["structure_sha256"], f"Source structure hash mismatch {identity}")
        checks = (("iptm", .7, True), ("peptide_plddt", 70., True), ("interface_pae", 10., False), ("minimum_heavyatom_distance", 1.5, True))
        failed = []
        for field, cutoff, lower in checks:
            value = float(row[field])
            require(math.isfinite(value), f"Nonfinite source quality {identity}")
            if (value < cutoff if lower else value > cutoff):
                failed.append(field)
        if int(row["contact_pairs_6A"]) <= 0:
            failed.append("contact_pairs_6A")
        require(row["structure_pass"].lower() == str(not failed).lower(), "Source quality flag disagrees with fixed gate")
        row.update(quality_pass=str(not failed), quality_fail_reasons=";".join(failed))
    require(identities == {(p, s, m) for p in REFERENCES for s in SEEDS for m in MODELS} and len(rows) == 60, "Incomplete 60-model panel")
    files = [source / name for name in ("reference_affinity.csv", "protocol.json", "analysis/model_results.csv")]
    files += [Path(__file__).resolve(), wrapper.resolve(), ROOT / "pipeline/verify_foldx_inputs.py", ROOT / "pipeline/verify_foldx_scores.py"] + [Path(r["complex_path"]) for r in rows]
    return rows, {str(p): sha(p) for p in files}


def check_hashes(hashes: dict[str, str], root: Path | None = None) -> None:
    for name, expected in hashes.items():
        path = root / name if root else Path(name)
        require(path.is_file() and sha(path) == expected, f"Frozen file missing/changed: {path}")


def prepare(source: Path = SOURCE, output: Path = OUTPUT, wrapper: Path | None = None) -> dict:
    """Prepare once; all later invocations independently validate the frozen files."""
    source, output = source.resolve(), output.resolve()
    wrapper = wrapper or ROOT / "scripts/run_foldx_affinity.sh"
    rows, source_hashes = validated_sources(source, wrapper)
    if (output / "protocol.json").exists():
        protocol = load_json(output / "protocol.json")
        require(load_json(output / "execution_state.json")["protocol_sha256"] == sha(output / "protocol.json"), "Frozen protocol changed")
        require(protocol["source_hashes"] == source_hashes, "Frozen source/script dependencies changed; refusing overwrite/resume")
        check_hashes(protocol["prepared_hashes"], output)
        require(protocol["scientific_protocol"] == scientific_protocol(), "Frozen scientific protocol changed")
        return protocol
    harmless = {"README.md", "run.log", "independent_input_verification.json"}
    require(not output.exists() or all(p.name in harmless and p.is_file() for p in output.iterdir()), "Unfrozen output has unknown artifacts; refusing overwrite")
    output.mkdir(parents=True, exist_ok=True)
    manifest, prepared_hashes = [], {}
    for row in rows:
        pid, seed, model = row["peptide_id"], int(row["seed"]), int(row["model_number"])
        job_id = f"{pid}_seed{seed}_model{model}"
        job = output / "jobs" / job_id
        job.mkdir(parents=True)
        content, counts = heavy_pdb(Path(row["complex_path"]).read_text(encoding="ascii"), {"A": RECEPTOR, "B": row["sequence"]})
        artifacts = {"input.pdb": content, "repair.cfg": "command=RepairPDB\npdb=input.pdb\n" + SETTINGS,
                     "analyse.cfg": "command=AnalyseComplex\npdb=input_Repair.pdb\nanalyseComplexChains=A,B\n" + SETTINGS}
        for name, value in artifacts.items():
            (job / name).write_text(value, encoding="ascii", newline="\n")
            prepared_hashes[str((job / name).relative_to(output))] = sha(job / name)
        manifest.append({"job_id": job_id, **{k: row[k] for k in MANIFEST_FIELDS if k in row}, "source_path": row["complex_path"],
                         "source_sha256": row["structure_sha256"], "input_path": str(job / "input.pdb"), "input_sha256": sha(job / "input.pdb"),
                         "repair_config_sha256": sha(job / "repair.cfg"), "analyse_config_sha256": sha(job / "analyse.cfg"), **counts})
    (output / "input_manifest.csv").write_text(csv_text(manifest, MANIFEST_FIELDS), encoding="utf-8", newline="\n")
    save_json(output / "preparation_verification.json", {"inputs": "PASS", "models_checked": 60, "heavy_atom_coordinates": "unchanged",
              "real_scoring": "NOT_RUN", "affinity_validation": "NOT_RUN", "panel": "DEVELOPMENT/SCREEN; previously inspected"})
    for name in ("input_manifest.csv", "preparation_verification.json"):
        prepared_hashes[name] = sha(output / name)
    protocol = {"state": "frozen_before_any_foldx_score", "source_dir": str(source), "wrapper": str(wrapper.resolve()),
                "scientific_protocol": scientific_protocol(), "source_hashes": source_hashes, "prepared_hashes": prepared_hashes}
    save_json(output / "protocol.json", protocol)
    save_json(output / "execution_state.json", {"state": "awaiting_foldx", "reason": "Official licensed FoldX executable not yet supplied", "jobs": {}, "protocol_sha256": sha(output / "protocol.json")})
    return protocol


def scientific_protocol() -> dict:
    return {"panel_role": "DEVELOPMENT/SCREEN; previously inspected; not blind external validation", "temperature_K": 298,
            "receptor_sequence": RECEPTOR, "peptide_length": 12, "seeds": list(SEEDS), "models": list(MODELS),
            "assay_approximation": "298 K approximates room temperature; 0.15 M represents major NaCl contribution", "pH": 7.4,
            "ionStrength_M": .15, "water": "-IGNORE", "screen": False, "commands": ["RepairPDB", "AnalyseComplex"], "chains": ["A", "B"],
            "noNterm": False, "noCterm": False, "terminal_treatment": "Free charged N terminus and C-terminal OXT retained on both chains",
            "units": "kcal/mol", "primary_minimum_quality_models_per_seed": 2, "minimum_matched_models": 2,
            "quality_gate": {"iptm_min": .7, "peptide_plddt_min": 70., "interface_pae_max": 10., "minimum_heavyatom_distance_min": 1.5, "contact_pairs_6A_strict_min": 0},
            "pdb_b_factors": "AFM pLDDT; not experimental temperature factors or affinity observations",
            "repair_documented_defaults": {"bFactorLimit": 0, "repair_Interface": "ALL", "real_tool_help_and_behavior": "inspect on first licensed run"},
            "pair_gain": "median(matched parent energy - variant energy)", "correlation": "Spearman(energy, log10(Kd_nM)); positive is correct",
            "strong_direction_cases": list(STRONG), "descriptive_only": ["P53_K24A"],
            "gate_each_seed": {"overall_rho_min": .6, "family_rho_strict_min": 0, "minimum_correct_directions": 6,
                               "required_correct": ["PMI_N8A", "P53_P27A"], "all_seven_comparable": True, "all_ten_peptides_eligible": True,
                               "global_requires_all_60_jobs_finished": True},
            "subsequent_steps": "Independent unused holdout and wet-lab Kd remain conditional; neither performed here"}


def parse_interaction(path: Path, expected_pdb: str | None = None) -> float:
    """Parse named official tabular fields; never infer energies from numeric log lines."""
    require(path.name.startswith(("Interaction_", "Summary_")) and path.suffix == ".fxout", "Not an AnalyseComplex interaction/summary table")
    headers, matches = None, []
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        cells = [cell.strip() for cell in line.split("\t")]
        if "Interaction Energy" in cells:
            require(headers is None and len(cells) == len(set(cells)), "Ambiguous repeated/duplicate header")
            require("Group1" in cells and "Group2" in cells, "Missing chain headers")
            headers = cells
            continue
        if headers is None or not line.strip():
            continue
        require(len(cells) == len(headers), "Malformed row after energy header")
        row = dict(zip(headers, cells))
        require({row["Group1"], row["Group2"]} == {"A", "B"}, "Wrong interaction chain pair")
        if expected_pdb is not None:
            require("Pdb" in row and Path(row["Pdb"]).stem == Path(expected_pdb).stem, "Wrong scored PDB")
        energy = float(row["Interaction Energy"])
        require(math.isfinite(energy), "Nonfinite interaction energy")
        matches.append(energy)
    require(headers is not None and len(matches) == 1, "Expected exactly one A/B Interaction Energy row")
    return matches[0]


def raw_energy(job: Path) -> tuple[float, list[str]]:
    interactions, summaries = sorted(job.glob("Interaction_*.fxout")), sorted(job.glob("Summary_*.fxout"))
    selected = interactions or summaries
    require(len(selected) == 1, "Missing/ambiguous AnalyseComplex energy files")
    energy = parse_interaction(selected[0], "input_Repair.pdb")
    for summary in summaries if interactions else []:
        require(math.isclose(parse_interaction(summary, "input_Repair.pdb"), energy, rel_tol=1e-8, abs_tol=1e-6), "Summary/interaction energy disagreement")
    return energy, [str(p) for p in selected]


def command_log(command: list[str], cwd: Path, stage: str, timeout: float) -> dict:
    started = time.monotonic()
    metadata = {"command": command, "cwd": str(cwd), "started_utc": datetime.now(timezone.utc).isoformat()}
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
    metadata["seconds"] = time.monotonic() - started
    save_json(cwd / f"{stage}.command.json", metadata)
    return metadata


def tool_manifest(output: Path, binary: Path, resource_dir: Path | None, timeout: float) -> dict:
    require(binary.is_absolute() and binary.is_file() and os.access(binary, os.X_OK), "FoldX must be an existing executable at an absolute path")
    binary = binary.resolve()
    directory = (resource_dir or binary.parent).resolve()
    require(directory.is_dir(), "Resource directory missing")
    resources = {}
    candidates = directory.rglob("*") if resource_dir else directory.glob("*")
    for path in sorted(candidates):
        if path.is_file() and (path.name.lower().startswith("rotabase") or resource_dir is not None and path.suffix.lower() in (".dat", ".tbl", ".param", ".params")):
            relative = path.relative_to(directory).as_posix()
            relative.encode("ascii")
            resources[relative] = {"source": str(path), "sha256": sha(path)}
    fingerprint = {"executable": str(binary), "executable_sha256": sha(binary), "resource_directory": str(directory), "resources": resources}
    frozen = output / "tool_manifest.json"
    if frozen.exists():
        manifest = load_json(frozen)
        for relative, resource in manifest["fingerprint"]["resources"].items():
            if Path(resource["source"]).parent == output / "tool_probe":
                require(sha(Path(resource["source"])) == resource["sha256"], "Generated help resource changed")
                fingerprint["resources"][relative] = resource
        require(manifest["fingerprint"] == fingerprint, "FoldX executable/resource fingerprint changed; refusing resume")
        check_hashes(manifest["capture_hashes"], output)
        return manifest
    probe = output / "tool_probe"
    probe.mkdir(exist_ok=True)
    log = command_log([str(binary), "--help"], probe, "help", min(timeout, 60))
    require(log["returncode"] == 0, f"FoldX --help failed: {log['error'] or log['returncode']}")
    help_text = (probe / "help.stdout.txt").read_text() + (probe / "help.stderr.txt").read_text()
    require(re.search(r"foldx", help_text, re.I) is not None, "Help does not identify FoldX")
    version4 = re.search(r"foldx\s*(?:suite\s*)?(?:version\s*[:=]?\s*)?4(?:\.|\b)", help_text, re.I) is not None
    generated = probe / "rotabase.txt"
    if generated.exists() and "rotabase.txt" not in resources:
        resources["rotabase.txt"] = {"source": str(generated), "sha256": sha(generated)}
        fingerprint["resources"] = resources
    require(not version4 or any(Path(name).name == "rotabase.txt" for name in resources), "FoldX v4 requires official rotabase.txt in --resource-dir")
    manifest = {"fingerprint": fingerprint, "help_command": log["command"], "help_identifies_version4": version4,
                "resource_policy": "Frozen supplied resources; v5 autogenerated parameters are bound to executable and hashed per job",
                "capture_hashes": {str(p.relative_to(output)): sha(p) for p in probe.iterdir() if p.is_file()}}
    save_json(frozen, manifest)
    return manifest


def job_signature(protocol: dict, tool: dict) -> str:
    return digest(json_text({"source": protocol["source_hashes"], "prepared": protocol["prepared_hashes"], "tool": tool}))


def job_artifacts(job: Path) -> dict:
    return {p.relative_to(job).as_posix(): sha(p) for p in sorted(job.rglob("*")) if p.is_file() and p.relative_to(job).as_posix() not in ("input.pdb", "repair.cfg", "analyse.cfg", "result.json")}


def validate_result(job: Path, row: dict, signature: str, expected_hash: str) -> dict:
    require(sha(job / "result.json") == expected_hash, f"Result record changed: {job.name}")
    result = load_json(job / "result.json")
    require(result["signature"] == signature and result["artifact_hashes"] == job_artifacts(job), f"Cache/raw artifacts changed: {job.name}")
    if result["status"] == "ok":
        heavy_pdb((job / "input_Repair.pdb").read_text(encoding="ascii"), {"A": RECEPTOR, "B": row["sequence"]})
        energy, _ = raw_energy(job)
        require(result["units"] == "kcal/mol" and result["interaction_energy_kcal_mol"] == energy, "Cached energy differs from raw output")
    else:
        require(result["status"] == "failed" and "interaction_energy_kcal_mol" not in result, "Failure cannot carry an energy")
    return result


def run(source: Path, output: Path, binary: Path | None, resource_dir: Path | None = None, timeout: float = 1800) -> dict:
    protocol = prepare(source, output)
    state = load_json(output / "execution_state.json")
    if binary is None or not binary.is_file():
        state.update(state="awaiting_foldx", reason=f"Official licensed FoldX executable unavailable: {binary}")
        save_json(output / "execution_state.json", state)
        raise ValueError(state["reason"] + "; no energies produced")
    try:
        tool = tool_manifest(output, binary, resource_dir, timeout)
    except (ValueError, OSError) as error:
        state.update(state="tool_initialization_failed", reason=str(error))
        save_json(output / "execution_state.json", state)
        raise
    signature = job_signature(protocol, tool)
    manifest = read_csv(output / "input_manifest.csv")
    priority = ["PMI_WT", "P53_WT", "PMI_N8A", "P53_P27A", "P53_L22A"]
    manifest.sort(key=lambda r: (priority.index(r["peptide_id"]) if r["peptide_id"] in priority else 99, r["peptide_id"], int(r["seed"]), int(r["model_number"])))
    state.update(state="running", reason="")
    save_json(output / "execution_state.json", state)
    for row in manifest:
        job = output / "jobs" / row["job_id"]
        if row["job_id"] in state["jobs"]:
            validate_result(job, row, signature, state["jobs"][row["job_id"]])
            continue
        require(not (job / "result.json").exists(), "Untracked result; refusing overwrite")
        # Preserve interrupted attempt output before trying again.
        old_artifacts = {p.name: sha(p) for p in job.iterdir() if p.is_file() and p.name not in ("input.pdb", "repair.cfg", "analyse.cfg", "result.json")}
        if old_artifacts:
            archive = job / "attempts" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%f")
            archive.mkdir(parents=True)
            for name in old_artifacts:
                shutil.move(str(job / name), str(archive / name))
        result = {"status": "failed", "signature": signature, "failure_reason": "", "units": "kcal/mol"}
        try:
            check_hashes(protocol["source_hashes"])
            check_hashes(protocol["prepared_hashes"], output)
            require(sha(Path(tool["fingerprint"]["executable"])) == tool["fingerprint"]["executable_sha256"], "Executable changed during run")
            for relative, resource in tool["fingerprint"]["resources"].items():
                require(sha(Path(resource["source"])) == resource["sha256"], "Tool resource changed during run")
                target = job / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(resource["source"], target)
            for stage in ("repair", "analyse"):
                log = command_log([str(binary.resolve()), "-f", f"{stage}.cfg"], job, stage, timeout)
                require(log["returncode"] == 0, f"{stage} failed: {log['error'] or log['returncode']}")
                for relative, resource in tool["fingerprint"]["resources"].items():
                    require(sha(Path(resource["source"])) == resource["sha256"] and sha(job / relative) == resource["sha256"], "Frozen tool resource changed during command")
                if stage == "repair":
                    heavy_pdb((job / "input_Repair.pdb").read_text(encoding="ascii"), {"A": RECEPTOR, "B": row["sequence"]})
            energy, paths = raw_energy(job)
            check_hashes(protocol["source_hashes"])
            check_hashes(protocol["prepared_hashes"], output)
            require(sha(Path(tool["fingerprint"]["executable"])) == tool["fingerprint"]["executable_sha256"], "Executable changed during command")
            result.update(status="ok", interaction_energy_kcal_mol=energy, raw_energy_paths=paths, repaired_path=str(job / "input_Repair.pdb"))
        except (ValueError, OSError) as error:
            result["failure_reason"] = str(error)
        result["artifact_hashes"] = job_artifacts(job)
        save_json(job / "result.json", result)
        state["jobs"][row["job_id"]] = sha(job / "result.json")
        save_json(output / "execution_state.json", state)
        print(f"{row['job_id']}: {result['status']}", flush=True)
    results = [load_json(output / "jobs" / row["job_id"] / "result.json") for row in manifest]
    state.update(state="complete" if all(r["status"] == "ok" for r in results) else "complete_with_failures",
                 scored_models=sum(r["status"] == "ok" for r in results), failed_models=sum(r["status"] == "failed" for r in results))
    save_json(output / "execution_state.json", state)
    return state


def paired_median_gain(parent: dict[int, float], variant: dict[int, float], minimum: int = 2) -> tuple[float | None, list[int]]:
    matched = sorted(parent.keys() & variant.keys())
    return (statistics.median(parent[m] - variant[m] for m in matched) if len(matched) >= minimum else None), matched


def analysis_tables(output: Path, protocol: dict) -> tuple[dict[str, list[dict]], dict]:
    manifest = read_csv(output / "input_manifest.csv")
    state, tool_path = load_json(output / "execution_state.json"), output / "tool_manifest.json"
    tool = load_json(tool_path) if tool_path.exists() else None
    signature = job_signature(protocol, tool) if tool else None
    model_rows, groups = [], {}
    for row in manifest:
        result = {"status": "not_run", "failure_reason": "No real FoldX result", "interaction_energy_kcal_mol": ""}
        if row["job_id"] in state["jobs"]:
            require(tool is not None, "Results lack frozen tool manifest")
            result = validate_result(output / "jobs" / row["job_id"], row, signature, state["jobs"][row["job_id"]])
        record = {**row, **{key: result.get(key, "") for key in ("status", "failure_reason", "interaction_energy_kcal_mol", "repaired_path")},
                  "units": "kcal/mol", "raw_energy_paths": ";".join(result.get("raw_energy_paths", []))}
        model_rows.append(record)
        groups.setdefault((row["peptide_id"], int(row["seed"])), []).append(record)
    summaries, directions, correlations, gates = [], [], [], []
    for scope in ("primary_quality", "all_models"):
        values = {}
        for (pid, seed), records in groups.items():
            selected = {int(r["model_number"]): float(r["interaction_energy_kcal_mol"]) for r in records if r["status"] == "ok" and (scope == "all_models" or r["quality_pass"] == "True")}
            eligible = len(selected) >= (2 if scope == "primary_quality" else 1)
            values[(pid, seed)] = selected if eligible else {}
            energies = list(selected.values())
            summaries.append({"scope": scope, "peptide_id": pid, "family": records[0]["family"], "seed": seed, "kd_nM": records[0]["kd_nM"],
                              "n_scored_selected": len(selected), "eligible": eligible, "median_energy_kcal_mol": statistics.median(energies) if eligible else "",
                              "minimum_energy_kcal_mol": min(energies) if energies else "", "maximum_energy_kcal_mol": max(energies) if energies else "",
                              "model_range_kcal_mol": max(energies) - min(energies) if energies else "", "n_failed": sum(r["status"] == "failed" for r in records),
                              "n_not_run": sum(r["status"] == "not_run" for r in records), "n_quality_failed": sum(r["quality_pass"] != "True" for r in records), "units": "kcal/mol"})
        for seed in SEEDS:
            seed_directions = []
            for pid, (_, kd, _, parent, _) in REFERENCES.items():
                if pid == parent:
                    continue
                gain, matched = paired_median_gain(values[(parent, seed)], values[(pid, seed)], 2 if scope == "primary_quality" else 1)
                expected = 1 if kd < REFERENCES[parent][1] else -1
                entry = {"scope": scope, "peptide_id": pid, "parent_id": parent, "seed": seed, "primary_direction_case": pid in STRONG,
                         "n_matched_models": len(matched), "matched_models": ";".join(map(str, matched)), "gain_kcal_mol": gain if gain is not None else "",
                         "experimental_direction": "stronger" if expected == 1 else "weaker", "comparable": gain is not None,
                         "direction_correct": gain * expected > 0 if gain is not None else "", "units": "kcal/mol"}
                directions.append(entry)
                if pid in STRONG:
                    seed_directions.append(entry)
            rhos = {}
            for family in ("overall", "PMI", "p53"):
                subset = [r for r in summaries if r["scope"] == scope and r["seed"] == seed and r["eligible"] and (family == "overall" or r["family"] == family)]
                rho = float(spearmanr([r["median_energy_kcal_mol"] for r in subset], [math.log10(float(r["kd_nM"])) for r in subset]).statistic) if len(subset) >= 2 and len({r["median_energy_kcal_mol"] for r in subset}) > 1 else None
                rhos[family] = rho
                correlations.append({"scope": scope, "seed": seed, "family": family, "n_peptides": len(subset), "rho_energy_log10Kd": rho if rho is not None else ""})
            if scope == "primary_quality":
                comparable = all(r["comparable"] for r in seed_directions)
                correct = sum(r["direction_correct"] is True for r in seed_directions)
                required = all(r["direction_correct"] is True for r in seed_directions if r["peptide_id"] in ("PMI_N8A", "P53_P27A"))
                n_eligible = sum(bool(values[(pid, seed)]) for pid in REFERENCES)
                all_ten = n_eligible == 10
                passed = all_ten and comparable and rhos["overall"] is not None and rhos["overall"] >= .6 and all(rhos[f] is not None and rhos[f] > 0 for f in ("PMI", "p53")) and correct >= 6 and required
                gates.append({"seed": seed, "all_seven_comparable": comparable, "all_ten_eligible": all_ten, "n_comparable_peptides": n_eligible,
                              "correct_of_seven": correct, "N8A_and_P27A_correct": required, "rho": rhos, "gate_pass": passed,
                              "status": "PASS" if passed else "FAIL" if comparable and all_ten else "NOT_EVALUABLE"})
    counts = {status: sum(r["status"] == status for r in model_rows) for status in ("ok", "failed", "not_run")}
    scoring_complete = state["state"] in ("complete", "complete_with_failures") and counts["not_run"] == 0
    return {"model_results.csv": model_rows, "sequence_summary.csv": summaries, "mutation_directions.csv": directions, "rank_correlations.csv": correlations}, {
        "panel_role": "DEVELOPMENT/SCREEN; not blind external validation", "counts": counts, "seed_gates": gates, "scoring_complete": scoring_complete,
        "continuation_pass": scoring_complete and all(g["gate_pass"] for g in gates),
        "real_scoring": "NOT_RUN" if counts["ok"] + counts["failed"] == 0 else "PARTIAL" if counts["not_run"] else "COMPLETE_WITH_FAILURES" if counts["failed"] else "COMPLETE",
        "holdout_validation": "NOT_RUN; conditional", "wet_lab_Kd": "NOT_RUN; conditional"}


def analyze(source: Path, output: Path) -> dict:
    protocol = prepare(source, output)
    tables, report = analysis_tables(output, protocol)
    directory = output / "analysis"
    directory.mkdir(exist_ok=True)
    for name, rows in tables.items():
        (directory / name).write_text(csv_text(rows), encoding="utf-8", newline="\n")
    save_json(directory / "decision.json", report)
    text = ["# FoldX MDM2 DEVELOPMENT/SCREEN 面板", "", "本面板已用于查看方向和选方法，不能称盲法外部验证。模型是结构不确定性样本，不能当作生物学重复。", "",
            f"真实评分状态：{report['real_scoring']}；成功 {report['counts']['ok']}，失败 {report['counts']['failed']}，未运行 {report['counts']['not_run']}。", "",
            "能量来自实际 AnalyseComplex A/B Interaction Energy（kcal/mol）；越低越有利。主统计要求每肽每 seed ≥2 个质量合格模型，配对方向要求 ≥2 个同号模型。gain 为配对 parent-minus-variant 的中位数。", ""]
    for gate in report["seed_gates"]:
        text.append(f"seed {gate['seed']}：{gate['status']}；可比肽 {gate['n_comparable_peptides']}/10；七项全部可比={gate['all_seven_comparable']}；方向正确 {gate['correct_of_seven']}/7；N8A/P27A 均正确={gate['N8A_and_P27A_correct']}；rho={gate['rho']}。")
    text += ["", "固定门槛：全部60个任务完成后，每 seed 十条肽均合格可比，总体 rho≥0.6，PMI/p53 rho各>0，至少6/7方向正确且含N8A/P27A，七项全部可比。K24A方向仅描述性，仍计入十肽总体相关。全模型统计、各模型范围、失败及无结果详见CSV。", "",
             f"继续判据通过={report['continuation_pass']}。后续未用作选方法的 holdout 和实验 Kd 均未完成，仍为条件性后续步骤。", ""]
    (directory / "report_zh.md").write_text("\n".join(text), encoding="utf-8")
    return report


def verify(source: Path, output: Path) -> dict:
    """Revalidate sources, raw files and recomputed published statistics even when pending."""
    protocol = prepare(source, output)
    tool_path = output / "tool_manifest.json"
    if tool_path.exists():
        tool = load_json(tool_path)
        fingerprint = tool["fingerprint"]
        require(sha(Path(fingerprint["executable"])) == fingerprint["executable_sha256"], "Frozen executable changed")
        for resource in fingerprint["resources"].values():
            require(sha(Path(resource["source"])) == resource["sha256"], "Frozen resource changed")
        check_hashes(tool["capture_hashes"], output)
    tables, report = analysis_tables(output, protocol)
    published = (output / "analysis/decision.json").exists()
    if published:
        for name, rows in tables.items():
            require((output / "analysis" / name).read_text(encoding="utf-8") == csv_text(rows), f"Published analysis differs from recomputed raw data: {name}")
        require(load_json(output / "analysis/decision.json") == report, "Published decision differs from raw-recomputed gate")
    result = {"inputs": "PASS", "frozen_hashes": "PASS", "raw_result_reparse": "PASS" if report["counts"]["ok"] else "NOT_APPLICABLE_NO_SCORED_MODELS",
              "published_tables": "PASS" if published else "NOT_PUBLISHED", **report}
    save_json(output / "verification.json", result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("prepare", "run", "analyze", "verify"))
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--foldx", type=Path)
    parser.add_argument("--resource-dir", type=Path)
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args()
    try:
        if args.command == "run":
            result = run(args.source.resolve(), args.output.resolve(), args.foldx, args.resource_dir, args.timeout)
        else:
            result = globals()[args.command](args.source.resolve(), args.output.resolve())
        print(json_text(result if args.command != "prepare" else {"inputs": "PASS", "models": 60, "real_scoring": "NOT_RUN_OR_SEE_EXECUTION_STATE"}))
    except (ValueError, OSError, KeyError) as error:
        parser.exit(2, f"FoldX benchmark stopped: {error}\n")


if __name__ == "__main__":
    main()
