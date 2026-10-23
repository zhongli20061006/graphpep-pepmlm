"""Excluded full-settings engine probe; it produces no affinity scores."""
import json
from pathlib import Path
import shutil
import sys
import time

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.backbone_relax import hydrate, relax


def main():
    output = OUT / "runtime_preflight"
    output.mkdir(exist_ok=True)
    source = ROOT / "output/foldx_matched_diagnostic_20261006/jobs/P53_P27A_seed20261011_model1/input_1_0.pdb"
    target = output / "input.pdb"
    if target.exists():
        assert target.read_bytes() == source.read_bytes()
    else:
        shutil.copyfile(source, target)
    started = time.monotonic()
    hydration = hydrate(target, output / "common", "excluded_runtime_preflight_P27A", "ETFSDLWKLLAE")
    reports = {}
    for condition in ("fixed", "flexible"):
        tick = time.monotonic()
        reports[condition] = relax(output / "common", output / condition, condition, "ETFSDLWKLLAE")
        reports[condition]["seconds"] = time.monotonic() - tick
        print(json.dumps(dict(condition=condition, seconds=reports[condition]["seconds"],
                              iterations=reports[condition]["iterations"], converged=reports[condition]["converged"],
                              frozen_max_A=reports[condition]["frozen_max_A"],
                              peptide_max_A=reports[condition]["B_backbone_max_A"])), flush=True)
    assert reports["fixed"]["frozen_max_A"] <= 1e-6 and reports["flexible"]["frozen_max_A"] <= 1e-6
    assert reports["fixed"]["B_backbone_max_A"] <= 1e-6 and reports["flexible"]["B_backbone_max_A"] > .01
    assert reports["fixed"]["start_positions_sha256"] == reports["flexible"]["start_positions_sha256"]
    report = dict(status="passed", excluded_from_formal_statistics=True, has_foldx_scores=False,
                  seconds=time.monotonic() - started, hydration=hydration, conditions=reports)
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
