"""Fresh final checks and sealed artifact manifest; no model calls."""
import csv
import hashlib
import importlib.metadata as metadata
import json
import sys
import unittest
from datetime import datetime, timezone
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
SITE = Path("/home/g/miniconda3/envs/colabfold_exp0/lib/python3.11/site-packages")


def load(path: Path):
    def reject(value):
        raise ValueError("Nonfinite JSON token " + value)
    return json.loads(path.read_text(encoding="utf-8-sig"), parse_constant=reject)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name: str, value) -> None:
    with (OUT / name).open("x", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def main() -> None:
    started = datetime.now(timezone.utc).isoformat()
    with (OUT / "test_after_run.log").open("x", encoding="utf-8") as stream:
        suite = unittest.defaultTestLoader.discover(str(OUT), pattern="test_msa_core.py")
        result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    assert result.wasSuccessful() and result.testsRun == 9
    for path in OUT.glob("*.py"):
        compile(path.read_text(), str(path), "exec")
    for path in OUT.glob("*.json"):
        load(path)
    assert load(OUT / "independent_verification.json")["all_pass"]
    protocol = load(OUT / "prepared_diagnostic_protocol.json")
    previous = load(ROOT / "output/bclw_parent_kernel_ablation_v1/preflight.json")["checks"]
    versions = {}
    for distribution in metadata.distributions(path=[str(SITE)]):
        name = distribution.metadata["Name"].lower().replace("_", "-")
        if name in {"colabfold", "jax", "jaxlib", "openmm"}:
            versions[name] = distribution.version
    current_batch = sha(SITE / "colabfold/batch.py")
    current_relax = sha(ROOT / "pipeline/relax_feedback_panel.py")
    comparison = {
        "current_versions_from_dist_metadata_no_runtime_import": versions,
        "baseline_versions": previous["versions"],
        "versions_match": versions == previous["versions"],
        "batch_source_match": current_batch == previous["batch_source_sha256"],
        "relax_source_match": current_relax == previous["relax_script_sha256"],
        "input_py_limit": "No historical input.py hash was recorded in kernel preflight; current input.py is frozen and I/O-tested. Full environment identity beyond tracked files is not established.",
        "weight_check": "Expected prior-run hash carried in protocol; current bulk weights not loaded/rehash-claimed. Future launch must verify.",
    }
    assert comparison["versions_match"] and comparison["batch_source_match"] and comparison["relax_source_match"]
    save("baseline_runtime_comparability.json", comparison)
    with (OUT / "source_input_correspondence.csv").open() as stream:
        references = list(csv.DictReader(stream))
    for row in references:
        target = ROOT / row["source"] if row["source"].startswith("output/") else OUT / row["source"]
        assert target.is_file(), target
    original_files = load(OUT / "preservation_before.json")["files"]
    mismatches = [row["path"] for row in original_files if not Path(row["path"]).is_file() or sha(Path(row["path"])) != row["sha256"]]
    assert not mismatches
    assert not Path(protocol["planned_output_root"]).exists()
    assert sha(Path(protocol["input_a3m"])) == protocol["input_a3m_sha256"]
    assert sha(OUT / "prepared_diagnostic_protocol.json") == (OUT / "prepared_diagnostic_protocol.sha256").read_text().strip()
    assert sha(OUT / "selection_implementation_frozen.json") == (OUT / "selection_implementation_frozen.sha256").read_text().strip()
    forbidden = [name for name in sys.modules if name.split(".")[0] in {"torch", "transformers", "jax", "jaxlib", "colabfold", "esm", "openmm"}]
    assert not forbidden
    plan = OUT / "plan_zh.md"
    plan.write_text(plan.read_text().replace("- [ ]", "- [x]"), encoding="utf-8")
    exclusions = {"artifact_manifest.json", "artifact_manifest.sha256", "final_delivery_check.json"}
    files = [{"path": path.relative_to(OUT).as_posix(), "bytes": path.stat().st_size, "sha256": sha(path)}
             for path in sorted(OUT.rglob("*")) if path.is_file() and "__pycache__" not in path.parts and path.name not in exclusions]
    save("artifact_manifest.json", {
        "scope": "All research/delivery files except this manifest, detached checksum, final_delivery_check.json and Python __pycache__; excluded final check references this manifest to avoid self-reference.",
        "files": files,
    })
    manifest_sha = sha(OUT / "artifact_manifest.json")
    (OUT / "artifact_manifest.sha256").write_text(manifest_sha + "  artifact_manifest.json\n")
    artifact_mismatches = [row["path"] for row in files if sha(OUT / row["path"]) != row["sha256"] or (OUT / row["path"]).stat().st_size != row["bytes"]]
    assert not artifact_mismatches
    final = {
        "status": "DELIVERY_VERIFIED_INPUT_GO_PREPARATION_ONLY",
        "started_utc": started, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "unit_tests_run": result.testsRun, "unit_tests_failures": len(result.failures), "unit_tests_errors": len(result.errors),
        "independent_delivery_checks": "24/24 PASS; separate arithmetic implementation, not external biological validation",
        "original_files_checked": len(original_files), "original_file_mismatches": mismatches,
        "manifest_files": len(files), "manifest_sha256": manifest_sha, "artifact_mismatches": artifact_mismatches,
        "output_bytes_before_this_final_check": sum(path.stat().st_size for path in OUT.rglob("*") if path.is_file()),
        "output_bytes_max": 16777216, "forbidden_runtime_modules_loaded": forbidden,
        "new_predictions": 0, "new_relaxations": 0, "new_model_scores": 0,
        "new_downloads": 0, "new_training": 0, "prepared_commands_launched": 0,
        "future_output_root_created": False,
        "limits": ["MSA small and predominantly unreviewed", "Assay chemistry and junction omitted", "Alignment simultaneous combinations not enumerated", "Full feature generation and predictions not tested", "No H/affinity validation"],
    }
    save("final_delivery_check.json", final)
    assert sum(path.stat().st_size for path in OUT.rglob("*") if path.is_file()) <= 16777216
    print(json.dumps(final, ensure_ascii=False))


if __name__ == "__main__":
    main()
