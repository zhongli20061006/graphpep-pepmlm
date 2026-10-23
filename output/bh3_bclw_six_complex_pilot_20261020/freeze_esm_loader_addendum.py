#!/usr/bin/env python3
"""Freeze the ESM loader compatibility amendment before any successful ESM score."""
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path


PILOT = Path(__file__).resolve().parent
ROOT = PILOT.parents[1]


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    config = json.loads((PILOT / "run_config_frozen_v1.json").read_text(encoding="utf-8"))
    execution = json.loads((PILOT / "run_execution_record_v1.json").read_text(encoding="utf-8"))
    original = PILOT / "esm_pll.py"
    compat = PILOT / "esm_pll_pytorch27.py"
    test = PILOT / "test_esm_pll_pytorch27_compat.py"
    checkpoint = Path(config["candidate_scorers"]["ESM2_sequence_prior"]["checkpoint"])
    failures = [e for e in execution["events"] if e["label"] == "esm2_sequence_pll" and e["returncode"] == 1]
    if sha(original) != config["candidate_scorers"]["ESM2_sequence_prior"]["code_sha256"]:
        raise SystemExit("original frozen ESM wrapper changed")
    if sha(checkpoint) != config["candidate_scorers"]["ESM2_sequence_prior"]["checkpoint_sha256"]:
        raise SystemExit("frozen ESM checkpoint changed")
    if not failures or (PILOT / "scores/esm_pll.csv").exists():
        raise SystemExit("expected failed initial load and no score file before compatibility retry")
    import torch
    unsafe = torch.serialization.get_unsafe_globals_in_checkpoint(str(checkpoint))
    if unsafe != ["argparse.Namespace"]:
        raise SystemExit(f"unexpected unsafe-global list; refuse allowlisting: {unsafe}")
    doc = {
        "schema": "bclw_six_complex_esm_loader_compat_addendum_v1",
        "status": "FROZEN_PRE_SCORE_NOT_RUN",
        "frozen_at_local": datetime.now().astimezone().isoformat(timespec="seconds"),
        "base_config_sha256": sha(PILOT / "run_config_frozen_v1.json"),
        "original_frozen_script_sha256": sha(original),
        "compatibility_script": "esm_pll_pytorch27.py",
        "compatibility_script_sha256": sha(compat),
        "compatibility_test": "test_esm_pll_pytorch27_compat.py",
        "compatibility_test_sha256": sha(test),
        "checkpoint_path": str(checkpoint),
        "checkpoint_sha256": sha(checkpoint),
        "root_cause": "fair-esm 2.0.0 load_model_and_alphabet_local calls torch.load without an explicit weights_only argument; PyTorch 2.7.1 defaults to weights_only=True and rejects argparse.Namespace stored in this local ESM checkpoint.",
        "static_unsafe_globals": unsafe,
        "compatibility_change": "Wrap only the existing local fair-esm loader in torch.serialization.safe_globals([argparse.Namespace]); preserve the default weights_only safe unpickler. No weights_only=False, external code, new checkpoint, PLL formula, sequence, or model setting.",
        "original_failed_attempt": {"label": failures[-1]["label"], "returncode": failures[-1]["returncode"],
                                    "log": failures[-1]["log"], "no_score_file_written": True},
        "tests_before_retry": {"test_count": 6, "passed": 6, "failed": 0,
                               "includes_allowlist_scope_and_context_test": True},
        "frozen_before_compatibility_retry": True,
        "retry_execution_status": "NOT_STARTED",
    }
    target = PILOT / "esm_loader_compat_addendum_v1.json"
    target.write_text(json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"addendum": str(target), "script_sha256": doc["compatibility_script_sha256"],
                      "test_sha256": doc["compatibility_test_sha256"], "checkpoint_sha256": doc["checkpoint_sha256"],
                      "unsafe_globals": unsafe}, ensure_ascii=True, indent=2))


if __name__ == "__main__":
    main()
