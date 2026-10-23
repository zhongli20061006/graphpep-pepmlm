"""Standalone figure from audited paired backbone summaries, with range labels."""
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path(__file__).resolve().parent
ANALYSIS = OUT / "analysis"


def main():
    audit = json.loads((ANALYSIS / "independent_audit.json").read_text())
    assert audit["status"] == "passed"
    source = ANALYSIS / "backbone_summary.csv"
    with source.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 24 and all(row["complete"] == "True" for row in rows)
    lookup = {(row["peptide_id"], row["condition"], int(row["seed"]), int(row["model_number"])): row for row in rows}
    panels = [("P53_P27A", "P27A (B11): expected stronger", "positive gain"),
              ("P53_L22A", "L22A (B6): expected weaker", "negative gain")]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.9), constrained_layout=True)
    identities = [(seed, model) for seed in (20261011, 20261012) for model in (1, 2, 3)]
    for ax, (pid, title, expected) in zip(axes, panels):
        for condition, color, offset in (("fixed", "#27639a", -.12), ("flexible", "#cf7527", .12)):
            medians, lows, highs = [], [], []
            for seed, model in identities:
                row = lookup[pid, condition, seed, model]
                medians.append(float(row["median_gain_kcal_mol"]))
                lows.append(float(row["minimum_gain_kcal_mol"]))
                highs.append(float(row["maximum_gain_kcal_mol"]))
            values = np.asarray(medians)
            ax.errorbar(np.arange(6) + offset, values,
                        yerr=[values - np.asarray(lows), np.asarray(highs) - values], fmt="o", color=color,
                        capsize=3, markersize=5, linewidth=1.5, label=condition.capitalize())
        ax.axhline(0, color="#555555", linestyle="--", linewidth=1)
        ax.axvline(2.5, color="#d4d4d4", linewidth=1)
        ax.set_xticks(np.arange(6), [f"{seed}\nm{model}" for seed, model in identities], fontsize=8)
        ax.set_title(title + "\n" + expected, fontsize=11)
        ax.set_ylabel("Paired WT - mutant Interaction Energy (kcal/mol)")
        ax.grid(axis="y", color="#ededed", linewidth=.8)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, loc="best")
    fig.suptitle("Fixed versus mobile peptide backbone\nDots: median of 5 paired repeats; bars: full range (no biological replication)", fontsize=12)
    for extension in ("png", "svg"):
        fig.savefig(ANALYSIS / ("fixed_flexible_comparison." + extension), dpi=200)
    plt.close(fig)
    metadata = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    audit_sha256=hashlib.sha256((ANALYSIS / "independent_audit.json").read_bytes()).hexdigest(),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    bars="minimum to maximum of five corresponding paired gains; not confidence intervals")
    (ANALYSIS / "figure_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
