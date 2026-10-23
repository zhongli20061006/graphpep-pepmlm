#!/usr/bin/env python3
"""Read-only launch gate for the frozen six-complex Bcl-w pilot."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


PILOT = Path(__file__).resolve().parent
ROOT = PILOT.parents[1]
CONFIG = PILOT / "run_config_frozen_v1.json"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    checks = []

    def check(label: str, path: Path, expected: str) -> None:
        found = path.is_file()
        actual = digest(path) if found else None
        checks.append({
            "label": label,
            "path": str(path),
            "exists": found,
            "expected_sha256": expected,
            "observed_sha256": actual,
            "pass": found and actual.lower() == expected.lower(),
        })

    for rel, expected in cfg["runtime_and_provenance"]["run_artifact_sha256"].items():
        path = ROOT / rel if rel.startswith("pipeline/") else PILOT / rel
        check(rel, path, expected)
    check(
        "frozen_input_mapping",
        ROOT / "output/bh3_bclw_source_audit_20261020/input_representation_frozen_pilot_v1.json",
        cfg["source_audit"]["frozen_input_mapping_sha256"],
    )
    check(
        "frozen_direction_criteria",
        ROOT / "output/reward_direction_benchmark_20261008/frozen_direction_criteria_v2.json",
        cfg["source_audit"]["direction_criteria_v2_sha256"],
    )
    check("AF2 model parameters", Path(cfg["prediction"]["model_parameter_file"]), cfg["prediction"]["model_parameter_sha256"])
    check("GraphPep checkpoint", Path(cfg["candidate_scorers"]["GraphPep"]["checkpoint"]), cfg["candidate_scorers"]["GraphPep"]["checkpoint_sha256"])
    check("ESM2 checkpoint", Path(cfg["candidate_scorers"]["ESM2_sequence_prior"]["checkpoint"]), cfg["candidate_scorers"]["ESM2_sequence_prior"]["checkpoint_sha256"])
    check("FoldX executable", Path(cfg["candidate_scorers"]["FoldX"]["binary"]), cfg["candidate_scorers"]["FoldX"]["binary_sha256"])
    for rel, expected in cfg["candidate_scorers"]["GraphPep"]["package_code_sha256"].items():
        path = PILOT / "graphpep_patches/resmap_compat.py" if rel == "compat_mapper" else Path("/home/g/GraphPep_v1.3") / rel
        check("GraphPep/" + rel, path, expected)

    output_paths = [PILOT / rel for rel in (
        "seed_20261008/structures", "seed_20261009/structures",
        "seed_20261008/relaxed", "seed_20261009/relaxed", "scores",
    )]
    output_state = {str(p.relative_to(PILOT)): p.exists() for p in output_paths}
    baseline_bytes = sum(p.stat().st_size for p in PILOT.rglob("*") if p.is_file())
    gpu_raw = subprocess.run(
        ["nvidia-smi", "--query-gpu=name,memory.total,memory.free,utilization.gpu", "--format=csv,noheader"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    gpu_apps = subprocess.run(
        ["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    active_gpu_apps = bool(gpu_apps) and "No running processes found" not in gpu_apps
    active = subprocess.run(["ps", "-eo", "pid=,args="], check=True, capture_output=True, text=True).stdout.splitlines()
    competing = [line.strip() for line in active if "colabfold_batch" in line and "pilot_runtime_preflight.py" not in line]
    disk_e = shutil.disk_usage("/mnt/e")
    disk_root = shutil.disk_usage("/")
    result = {
        "schema": "bclw_six_complex_pilot_launch_preflight_v1",
        "config_sha256": digest(CONFIG),
        "hash_checks": checks,
        "hashes_pass": all(c["pass"] for c in checks),
        "gpu_snapshot": gpu_raw,
        "active_gpu_compute_apps": gpu_apps,
        "active_gpu_compute_apps_found": active_gpu_apps,
        "other_colabfold_processes": competing,
        "disk_bytes_free": {"workspace_mount_E": disk_e.free, "wsl_root": disk_root.free},
        "pilot_output_baseline_bytes": baseline_bytes,
        "planned_output_dirs": output_state,
        "expected_output_dirs_clean": not any(output_state.values()),
        "minimum_E_headroom_pass": disk_e.free >= int(cfg["budgets"]["max_new_output_size_GiB"] * 1024**3),
        "execution_status_before_launch": "NOT_STARTED",
    }
    out = PILOT / "launch_preflight_v1.json"
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "preflight_file": str(out),
        "hash_checks": len(checks),
        "hashes_pass": result["hashes_pass"],
        "gpu_snapshot": gpu_raw,
        "active_gpu_compute_apps": gpu_apps,
        "other_colabfold_processes": competing,
        "free_E_GiB": round(disk_e.free / 1024**3, 2),
        "free_WSL_root_GiB": round(disk_root.free / 1024**3, 2),
        "baseline_bytes": baseline_bytes,
        "output_dirs_clean": result["expected_output_dirs_clean"],
        "minimum_E_headroom_pass": result["minimum_E_headroom_pass"],
    }, ensure_ascii=True, indent=2))
    return 0 if (result["hashes_pass"] and result["expected_output_dirs_clean"]
                 and result["minimum_E_headroom_pass"] and not active_gpu_apps and not competing) else 2


if __name__ == "__main__":
    raise SystemExit(main())
