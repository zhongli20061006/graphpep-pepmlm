#!/usr/bin/env python3
"""Run one frozen-pilot command with cumulative GPU-time and output-size caps."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


PILOT = Path(__file__).resolve().parent
CONFIG = PILOT / "run_config_frozen_v1.json"
PREFLIGHT = PILOT / "launch_preflight_v1.json"
RECORD = PILOT / "run_execution_record_v1.json"
GIB = 1024**3


def now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def byte_size() -> int:
    return sum(p.stat().st_size for p in PILOT.rglob("*") if p.is_file())


def save(data: dict) -> None:
    tmp = RECORD.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    tmp.replace(RECORD)


def terminate_group(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        proc.wait(timeout=20)
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", required=True)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--env", action="append", default=[], help="KEY=VALUE environment override")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command and args.command[0] == "--" else args.command
    if not command:
        parser.error("a command is required after --")

    cfg = json.loads(CONFIG.read_text(encoding="utf-8"))
    preflight = json.loads(PREFLIGHT.read_text(encoding="utf-8"))
    if not preflight.get("hashes_pass") or not preflight.get("expected_output_dirs_clean"):
        raise SystemExit("preflight did not pass; refusing to launch")
    if digest(CONFIG) != preflight["config_sha256"]:
        raise SystemExit("frozen config changed after preflight; refusing to launch")

    if RECORD.exists():
        record = json.loads(RECORD.read_text(encoding="utf-8"))
        if record.get("config_sha256") != digest(CONFIG):
            raise SystemExit("execution record belongs to a different frozen config")
    else:
        record = {
            "schema": "bclw_six_complex_pilot_execution_record_v1",
            "status": "RUNNING",
            "started_at": now(),
            "config_sha256": digest(CONFIG),
            "preflight_sha256": digest(PREFLIGHT),
            "output_baseline_bytes": preflight["pilot_output_baseline_bytes"],
            "gpu_wall_seconds": 0.0,
            "output_cap_bytes": int(cfg["budgets"]["max_new_output_size_GiB"] * GIB),
            "events": [],
        }
    elapsed_gpu = float(record["gpu_wall_seconds"])
    gpu_cap = float(cfg["budgets"]["max_gpu_hours"]) * 3600
    remaining = gpu_cap - elapsed_gpu
    if args.gpu and remaining <= 0:
        raise SystemExit("cumulative GPU wall-time budget is exhausted")

    env = os.environ.copy()
    for item in args.env:
        key, sep, value = item.partition("=")
        if not sep or not key:
            parser.error(f"invalid --env value: {item}")
        env[key] = value

    logs = PILOT / "logs"
    logs.mkdir(exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in args.label)
    logfile = logs / f"{safe}.log"
    event = {"label": args.label, "command": command, "gpu_metered": args.gpu, "started_at": now(), "log": str(logfile.relative_to(PILOT))}
    started = time.monotonic()
    reason = None
    rc = None
    with logfile.open("ab", buffering=0) as log:
        log.write((f"\n[{now()}] START {args.label}\nCMD: {command!r}\n").encode())
        proc = subprocess.Popen(command, cwd=PILOT.parents[1], env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        event["pid"] = proc.pid
        while proc.poll() is None:
            time.sleep(5)
            out_delta = byte_size() - int(record["output_baseline_bytes"])
            if out_delta > int(record["output_cap_bytes"]):
                reason = "output_size_limit_exceeded"
                terminate_group(proc)
                break
            if args.gpu and time.monotonic() - started >= remaining:
                reason = "gpu_wall_time_limit_reached"
                terminate_group(proc)
                break
            if int(time.monotonic() - started) and int(time.monotonic() - started) % 30 < 5:
                print(f"{args.label}: running; elapsed={time.monotonic()-started:.0f}s; output_delta={out_delta/GIB:.3f} GiB", flush=True)
        rc = proc.wait()
        log.write((f"[{now()}] END {args.label} rc={rc} reason={reason}\n").encode())

    duration = time.monotonic() - started
    if args.gpu:
        record["gpu_wall_seconds"] = elapsed_gpu + duration
    event.update({"ended_at": now(), "duration_seconds": round(duration, 3), "returncode": rc, "stop_reason": reason,
                  "output_bytes_delta": byte_size() - int(record["output_baseline_bytes"])})
    record["events"].append(event)
    record["status"] = "RUNNING" if rc == 0 and reason is None else "COMMAND_FAILED_OR_STOPPED"
    record["updated_at"] = now()
    save(record)
    print(json.dumps(event, ensure_ascii=False, indent=2))
    return 0 if rc == 0 and reason is None else (rc if rc not in (None, 0) else 124)


if __name__ == "__main__":
    raise SystemExit(main())
