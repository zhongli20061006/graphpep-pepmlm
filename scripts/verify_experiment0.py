#!/usr/bin/env python3
"""Verify local artifacts and the author example without cloud predictions."""
from __future__ import annotations
import ast
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd

WORKSPACE = Path("/mnt/e/Users/graphpep与pepmlm")
WSL_ROOT = Path("/home/g/pepmlm_pipeline")
PILOT = WORKSPACE / "output/experiment0"
GRAPHPEP = Path("/home/g/GraphPep_v1.3")
sys.path.insert(0, str(WORKSPACE / "pipeline"))
import step3_graphpep as gp


def main() -> int:
    config = json.loads((PILOT / "pilot_config.json").read_text())
    assert config["canonical_sequence_verified"]
    assert config["target_length"] == 858 and config["receptor_length"] == 699
    assert config["receptor_range_1based_inclusive"] == [20, 718]
    notebook = json.loads((PILOT / "Experiment0_ColabFold.ipynb").read_text())
    checks = []
    for cell in notebook["cells"]:
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        if "%%bash" in source:
            shell = source.split("%%bash", 1)[1].split("\n", 1)[1]
            result = subprocess.run(["bash", "-n"], input=shell, text=True, capture_output=True)
            assert result.returncode == 0, result.stderr
        else:
            ast.parse(source)
        checks.append(cell["metadata"]["id"])
    fasta = (PILOT / "pilot_complexes.fasta").read_text().splitlines()
    assert len(fasta) == 4
    for index, row in enumerate(config["candidates"]):
        rec, pep = fasta[index * 2 + 1].split(":")
        assert rec == config["receptor_sequence"] and pep == row["peptide_sequence"]

    # Synchronize only changed pipeline files, preserving prior versions separately.
    backup = WSL_ROOT / "work/experiment0_backups_2026-10-02"
    backup.mkdir(parents=True, exist_ok=True)
    for name in ["step3_graphpep.py", "step4_analyze.py", "experiment0.py", "colabfold_batch_reference.ipynb"]:
        source, target = WORKSPACE / "pipeline" / name, WSL_ROOT / "pipeline" / name
        if target.exists() and not (backup / name).exists():
            shutil.copy2(target, backup / name)
        target.write_text(source.read_text())
        if name.endswith(".py"):
            ast.parse(target.read_text())

    source_files = [GRAPHPEP / "GraphPep.sh", GRAPHPEP / "bin/core.py", GRAPHPEP / "bin/data.py"]
    before = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files}
    control = PILOT / "software_control"
    control.mkdir(exist_ok=True)
    example = GRAPHPEP / "example/1B9J"
    rec_models = gp.read_models(example / "protein.pdb")
    pep_models = gp.read_models(example / "decoys.pdb")
    assert len(rec_models) == 1 and pep_models
    peptide_seq = gp.split_chain(pep_models[0], "B")[1]
    (control / "CONTROL_1B9J_complex.pdb").write_text("".join(rec_models[0] + pep_models[0]) + "END\n")
    manifest = pd.DataFrame([{"peptide_id": "CONTROL_1B9J", "peptide_sequence": peptide_seq,
                              "source": "software_control", "model_id": "author_decoy_1"}])
    pair = gp.build_inputs(manifest, str(control))[0]
    result = gp.score_pair(pair, GRAPHPEP, WSL_ROOT / "pipeline/graphpep_patches/resmap_compat.py", .7)
    reference = pd.read_csv(example / "score.csv")
    expected = float(reference.iloc[0]["PRED_SCORE"])
    assert abs(result["binding_score"] - expected) <= 5e-5, (result["binding_score"], expected)
    after = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files}
    assert before == after, "GraphPep installation source was changed"
    pd.DataFrame([result]).to_csv(control / "control_score.csv", index=False)
    verification = {"local_artifacts_verified": True, "cloud_prediction_started": False,
                    "notebook_syntax_checked_cells": checks, "canonical_target_verified": True,
                    "n_peptides": 2, "complex_length": 709,
                    "software_control_score": result["binding_score"], "author_reference_score": expected,
                    "control_n_contacts": result["n_contacts"], "graphpep_source_unchanged": before == after,
                    "graphpep_source_sha256": before,
                    "note": "Software control is not a generated candidate or a cloud prediction."}
    (PILOT / "local_verification.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(json.dumps(verification, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
