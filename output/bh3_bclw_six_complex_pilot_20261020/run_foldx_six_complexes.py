#!/usr/bin/env python3
"""Run the frozen RepairPDB/AnalyseComplex protocol on six relaxed complexes."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path


PILOT = Path(__file__).resolve().parent
ROOT = PILOT.parents[1]
CONFIG = json.loads((PILOT / "run_config_frozen_v1.json").read_text(encoding="utf-8"))
TEMPLATE = json.loads((PILOT / "foldx_protocol_template.json").read_text(encoding="utf-8"))
FOLDX = Path(TEMPLATE["foldx_executable"])


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_interaction(path: Path, expected_pdb: str) -> float:
    headers, values = None, []
    for line in path.read_text(encoding="utf-8", errors="strict").splitlines():
        cells = [cell.strip() for cell in line.split("\t")]
        if "Interaction Energy" in cells:
            if headers is not None or len(cells) != len(set(cells)):
                raise ValueError("ambiguous or duplicate Interaction Energy header")
            headers = cells
            if not {"Pdb", "Group1", "Group2", "Interaction Energy"}.issubset(headers):
                raise ValueError("required FoldX interaction columns are missing")
            continue
        if headers is None or not line.strip():
            continue
        if len(cells) != len(headers):
            raise ValueError("malformed FoldX interaction row")
        row = dict(zip(headers, cells))
        if {row["Group1"], row["Group2"]} != {"A", "B"}:
            raise ValueError("FoldX did not report the A/B interaction")
        if Path(row["Pdb"]).stem != Path(expected_pdb).stem:
            raise ValueError("FoldX interaction output names an unexpected PDB")
        value = float(row["Interaction Energy"])
        if not math.isfinite(value):
            raise ValueError("nonfinite Interaction Energy")
        values.append(value)
    if headers is None or len(values) != 1:
        raise ValueError("expected exactly one A/B Interaction Energy row")
    return values[0]


def run_stage(job: Path, stage: str) -> dict:
    command = [str(FOLDX), "-f", f"{stage}.cfg"]
    try:
        result = subprocess.run(command, cwd=job, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", timeout=900,
                                env={**__import__("os").environ, "LC_ALL": "C", "LANG": "C"}, shell=False)
        (job / f"{stage}.stdout.txt").write_text(result.stdout, encoding="utf-8")
        (job / f"{stage}.stderr.txt").write_text(result.stderr, encoding="utf-8")
        meta = {"command": command, "returncode": result.returncode, "error": ""}
    except (OSError, subprocess.TimeoutExpired) as exc:
        stdout = getattr(exc, "stdout", "") or ""
        stderr = getattr(exc, "stderr", "") or ""
        if isinstance(stdout, bytes): stdout = stdout.decode(errors="replace")
        if isinstance(stderr, bytes): stderr = stderr.decode(errors="replace")
        (job / f"{stage}.stdout.txt").write_text(stdout, encoding="utf-8")
        (job / f"{stage}.stderr.txt").write_text(stderr, encoding="utf-8")
        meta = {"command": command, "returncode": None, "error": str(exc)}
    (job / f"{stage}.command.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    if meta["returncode"] != 0:
        raise RuntimeError(f"{stage} failed: {meta['error'] or meta['returncode']}")
    return meta


def main() -> int:
    if not FOLDX.is_file() or sha(FOLDX) != CONFIG["candidate_scorers"]["FoldX"]["binary_sha256"]:
        raise SystemExit("FoldX executable is missing or its frozen hash changed")
    rows = []
    with (PILOT / "graphpep_inputs.csv").open(encoding="utf-8-sig", newline="") as f:
        manifest = list(csv.DictReader(f))
    for row in manifest:
        ident, model_id = row["peptide_id"], row["model_id"]
        seed = model_id.removeprefix("seed")
        job = PILOT / "foldx_jobs" / f"{ident}_{model_id}"
        job.mkdir(parents=True, exist_ok=True)
        source = Path(row["complex_path"])
        result = {"peptide_id": ident, "model_id": model_id, "seed": int(seed), "source_pdb": str(source),
                  "source_sha256": sha(source), "status": "failed", "failure_reason": ""}
        try:
            if any(job.iterdir()):
                raise RuntimeError("job directory is not empty; refusing overwrite")
            shutil.copyfile(source, job / "input.pdb")
            (job / "repair.cfg").write_text(TEMPLATE["repair_config_template"], encoding="ascii", newline="\n")
            (job / "analyse.cfg").write_text(TEMPLATE["analyse_config_template"], encoding="ascii", newline="\n")
            # The authorized install is checked for optional Rotabase resources; this install contains only the binary.
            for resource in FOLDX.parent.glob("rotabase*"):
                if resource.is_file():
                    shutil.copyfile(resource, job / resource.name)
            result["input_sha256"] = sha(job / "input.pdb")
            result["repair_config_sha256"] = sha(job / "repair.cfg")
            result["analyse_config_sha256"] = sha(job / "analyse.cfg")
            run_stage(job, "repair")
            repaired = job / "input_Repair.pdb"
            if not repaired.is_file():
                raise RuntimeError("RepairPDB did not produce input_Repair.pdb")
            run_stage(job, "analyse")
            interactions = sorted(job.glob("Interaction_*.fxout"))
            if len(interactions) != 1:
                raise RuntimeError(f"expected one Interaction_*.fxout; got {len(interactions)}")
            energy = parse_interaction(interactions[0], "input_Repair.pdb")
            summaries = sorted(job.glob("Summary_*.fxout"))
            if summaries:
                if len(summaries) != 1 or not math.isclose(parse_interaction(summaries[0], "input_Repair.pdb"), energy, rel_tol=1e-8, abs_tol=1e-6):
                    raise RuntimeError("FoldX Summary and Interaction outputs disagree")
            result.update(status="ok", interaction_energy_kcal_mol=energy,
                          interaction_file=str(interactions[0].relative_to(PILOT)),
                          repaired_pdb_sha256=sha(repaired))
        except Exception as exc:  # each planned cell remains represented even on tool failure
            result["failure_reason"] = f"{type(exc).__name__}: {exc}"
        result["job_directory"] = str(job.relative_to(PILOT))
        (job / "result.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        rows.append(result)
        print(f"{ident} {model_id}: {result['status']} {result.get('interaction_energy_kcal_mol', result['failure_reason'])}", flush=True)

    out = PILOT / "scores/foldx.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = sorted({key for row in rows for key in row})
    with out.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    summary = {"schema": "bclw_six_complex_foldx_execution_v1", "expected_rows": 6,
               "observed_rows": len(rows), "ok_rows": sum(r["status"] == "ok" for r in rows),
               "failed_rows": sum(r["status"] != "ok" for r in rows),
               "interpretation": "diagnostic proxy only; all six structures failed the frozen AFM confidence gate",
               "rows": rows}
    (PILOT / "scores/foldx_execution.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0 if len(rows) == 6 else 2


if __name__ == "__main__":
    raise SystemExit(main())
