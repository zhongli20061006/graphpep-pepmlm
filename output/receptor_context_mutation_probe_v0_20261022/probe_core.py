"""Standard-library-only core for the frozen receptor-context probe.

This module performs bookkeeping and arithmetic only; it never loads a model.
"""
from __future__ import annotations

from collections import OrderedDict
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import stat
import subprocess
import time
from typing import Any, Mapping, Sequence


EPSILON = 1e-4
SITE_COUNT = 15
JOB_CAP = 60


def sha256(path: str | os.PathLike[str]) -> str:
    """Return the SHA-256 digest of a file without reading it all into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: str | os.PathLike[str]) -> Any:
    with Path(path).open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def write_json(path: str | os.PathLike[str], obj: Any) -> None:
    """Write reproducible, readable JSON, rejecting nonfinite numeric values."""
    payload = json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)
    Path(path).write_text(payload + "\n", encoding="utf-8")


def directory_bytes(root: str | os.PathLike[str]) -> int:
    """Sum regular-file bytes below root, without following symbolic links."""
    total = 0
    root = Path(root)
    if not root.is_dir():
        raise ValueError("output root must be an existing directory")
    for directory, _, filenames in os.walk(root, followlinks=False):
        for filename in filenames:
            try:
                info = (Path(directory) / filename).stat(follow_symlinks=False)
            except FileNotFoundError:
                # A worker may replace a temporary file during a scan.
                continue
            if stat.S_ISREG(info.st_mode):
                total += info.st_size
    return total


def mask_tokens(
    token_ids: Sequence[int],
    receptor_ids: Sequence[int],
    peptide_ids: Sequence[int],
    site: int,
    cls_id: int,
    eos_id: int,
    mask_id: int,
) -> tuple[list[int], int]:
    """Validate CLS+R+P+EOS and replace exactly one 1-based peptide site."""
    receptor = list(receptor_ids)
    peptide = list(peptide_ids)
    expected = [cls_id, *receptor, *peptide, eos_id]
    if list(token_ids) != expected:
        raise ValueError("token_ids must exactly equal CLS + receptor + peptide + EOS")
    if isinstance(site, bool) or not isinstance(site, int) or not 1 <= site <= len(peptide):
        raise ValueError("site must be a 1-based position within the peptide")
    if mask_id in expected:
        raise ValueError("unmasked input must not already contain the mask token")
    index = len(receptor) + site
    masked = expected.copy()
    masked[index] = mask_id
    return masked, index


def contrast(
    target: Mapping[str, float],
    control_1: Mapping[str, float],
    control_2: Mapping[str, float],
    wt: str,
    alt: str,
) -> dict[str, float | bool]:
    """Compute conditional log-odds differences under frozen strict gates."""
    differences = []
    for context in (target, control_1, control_2):
        wt_logp, alt_logp = float(context[wt]), float(context[alt])
        if not math.isfinite(wt_logp) or not math.isfinite(alt_logp):
            raise ValueError("log probabilities must be finite")
        differences.append(alt_logp - wt_logp)
    d_target, d_control_1, d_control_2 = differences
    h_control_1 = d_target - d_control_1
    h_control_2 = d_target - d_control_2
    h = d_target - (0.5 * d_control_1 + 0.5 * d_control_2)
    if not all(math.isfinite(value) for value in (
        d_target, d_control_1, d_control_2, h_control_1, h_control_2, h
    )):
        raise ValueError("contrast arithmetic must remain finite")
    disagreement = (
        h_control_1 > EPSILON and h_control_2 < -EPSILON
    ) or (
        h_control_2 > EPSILON and h_control_1 < -EPSILON
    )
    return {
        "d_target": d_target,
        "d_control_1": d_control_1,
        "d_control_2": d_control_2,
        "h_control_1": h_control_1,
        "h_control_2": h_control_2,
        "h": h,
        "passes_filter": all(value > EPSILON for value in (
            d_target, h_control_1, h_control_2
        )),
        "control_sign_disagreement": disagreement,
    }


def select_suggestions(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    """Keep up to six distinct sites per parent, returning the original rows."""
    parents: OrderedDict[str, list[dict[str, Any]]] = OrderedDict()
    for row in rows:
        if row.get("excluded_historical", False) or not row.get("passes_filter", False):
            continue
        if not math.isfinite(float(row["h"])):
            raise ValueError("suggestion h values must be finite")
        parents.setdefault(row["parent_id"], []).append(row)
    selected = []
    for candidates in parents.values():
        sites = set()
        for row in sorted(candidates, key=lambda item: (
            -float(item["h"]), item["site"], item["alt"]
        )):
            if row["site"] in sites:
                continue
            selected.append(row)
            sites.add(row["site"])
            if len(sites) == 6:
                break
    return selected


def build_jobs(
    inputs: Sequence[Mapping[str, Any]], context_ids: Sequence[str]
) -> list[dict[str, Any]]:
    """Freeze 45 base jobs and 15 target repeats in reverse target order."""
    contexts = list(context_ids)
    if (len(contexts) != 3
            or any(not isinstance(context, str) or not context for context in contexts)
            or len(set(contexts)) != 3):
        raise ValueError("exactly three distinct nonempty context IDs are required")
    sites = []
    parent_ids = set()
    for item in inputs:
        parent_id = item["parent_id"]
        if not isinstance(parent_id, str) or not parent_id or parent_id in parent_ids:
            raise ValueError("parent IDs must be distinct nonempty strings")
        parent_ids.add(parent_id)
        eligible = list(item["eligible_positions"])
        if (any(isinstance(site, bool) or not isinstance(site, int) or site < 1
                for site in eligible) or len(set(eligible)) != len(eligible)):
            raise ValueError("eligible positions must be distinct positive integers")
        sites.extend((parent_id, site) for site in eligible)
    if len(sites) != SITE_COUNT:
        raise ValueError("this frozen experiment requires exactly 15 eligible sites")
    jobs = []
    for parent_id, site in sites:
        for context in contexts:
            jobs.append({
                "job_id": f"base_{len(jobs) + 1:03d}",
                "parent_id": parent_id,
                "site": site,
                "context_id": context,
                "phase": "base",
            })
    target_jobs = [job for job in jobs if job["context_id"] == contexts[0]]
    for index, base in enumerate(reversed(target_jobs), start=1):
        jobs.append({
            "job_id": f"repeat_{index:03d}",
            "parent_id": base["parent_id"],
            "site": base["site"],
            "context_id": base["context_id"],
            "phase": "repeat",
        })
    if len(jobs) != JOB_CAP:
        raise ValueError("the frozen job inventory must contain exactly 60 jobs")
    return jobs


def _kill_own_group(process: subprocess.Popen) -> None:
    """The session created by Popen has the worker PID as its process-group ID."""
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def run_bounded(
    command: Sequence[str],
    root: str | os.PathLike[str],
    seconds: float,
    max_bytes: int,
    log_path: str | os.PathLike[str],
) -> dict[str, Any]:
    """Run one Linux process group with wall-time and total-output-byte caps."""
    if os.name != "posix":
        raise RuntimeError("the bounded worker must be run on Linux")
    if not math.isfinite(seconds) or seconds <= 0:
        raise ValueError("seconds must be finite and positive")
    if isinstance(max_bytes, bool) or not isinstance(max_bytes, int) or max_bytes < 0:
        raise ValueError("max_bytes must be a nonnegative integer")
    if isinstance(command, (str, bytes)) or not command:
        raise ValueError("command must be a nonempty argument sequence")
    root = Path(root).resolve()
    if not root.is_dir():
        raise ValueError("output root must be an existing directory")
    log_path = Path(log_path).resolve()
    if not log_path.is_relative_to(root):
        raise ValueError("log_path must be within the output root")

    started = time.monotonic()
    process = None
    normal_exit = False
    status = None
    try:
        with log_path.open("wb") as log:
            process = subprocess.Popen(
                list(command), cwd=root, stdout=log, stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            while True:
                # Poll before scanning so even a rapidly exited worker has its
                # final output counted before completion is considered.
                returncode = process.poll()
                output_bytes = directory_bytes(root)
                elapsed = time.monotonic() - started
                if elapsed >= seconds:
                    status = "wall_cap"
                elif output_bytes > max_bytes:
                    status = "disk_cap"
                elif returncode is not None:
                    status = "completed" if returncode == 0 else "worker_failed"
                    break
                if status is not None:
                    _kill_own_group(process)
                    break
                time.sleep(0.05)
            process.wait()
        normal_exit = True
    finally:
        if process is not None:
            if not normal_exit or process.poll() is None:
                _kill_own_group(process)
            process.wait()

    output_bytes = directory_bytes(root)
    wall_seconds = time.monotonic() - started
    # Closing handles and collecting the final bytes also belong to the budget.
    # Recheck both caps even when the worker exited between two polling ticks.
    if status in ("completed", "worker_failed"):
        if wall_seconds >= seconds:
            status = "wall_cap"
        elif output_bytes > max_bytes:
            status = "disk_cap"
        if status in ("wall_cap", "disk_cap"):
            _kill_own_group(process)
            process.wait()
            output_bytes = directory_bytes(root)
            wall_seconds = time.monotonic() - started
    return {
        "status": status,
        "returncode": process.returncode,
        "wall_seconds": wall_seconds,
        "output_bytes": output_bytes,
        "pid": process.pid,
    }
