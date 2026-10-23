"""Export small-sample feedback comparisons as PNG and SVG."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/feedback_panel_mdm2_20261005"
ARMS = ["native_guided", "native_random", "predicted_guided", "predicted_random"]
LABELS = ["Native\nguided", "Native\nrandom", "Predicted\nguided", "Predicted\nrandom"]
PARENTS = ["P53WT_1YCR", "PMI_3EQS", "P53LIKE_3G03"]
COLORS = ["#087f8c", "#999999", "#c87028", "#999999"]


def main():
    seq = pd.read_csv(OUT / "analysis/sequence_summary.csv")
    arm = pd.read_csv(OUT / "analysis/arm_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.3), sharey=True, constrained_layout=True)
    for axis, parent in zip(axes, PARENTS):
        group = seq[seq.parent_id == parent]
        for x, name in enumerate(ARMS):
            points = group[group.arm == name]
            usable = points[points.median_gain_vs_predicted_parent.notna()]
            offsets = np.linspace(-.12, .12, len(usable)) if len(usable) > 1 else np.asarray([0.])
            if len(usable):
                axis.scatter(x + offsets, usable.median_gain_vs_predicted_parent, s=45,
                             color=COLORS[x], edgecolor="white", linewidth=.7, zorder=3)
                axis.plot([x-.22, x+.22], [usable.median_gain_vs_predicted_parent.median()]*2,
                          color=COLORS[x], linewidth=2, zorder=4)
            row = arm[(arm.parent_id == parent) & (arm.arm == name)].iloc[0]
            axis.text(x, .98, f"{int(row.children_improved)}/3", transform=axis.get_xaxis_transform(), ha="center", va="top", fontsize=10)
        axis.axhline(0, color="#444444", linewidth=.8)
        axis.axhline(.05, color="#444444", linewidth=.8, linestyle="--")
        axis.set_xticks(range(4), LABELS, fontsize=9)
        axis.set_title(parent, fontsize=11, pad=17)
        axis.set_xlim(-.5, 3.5)
        axis.grid(axis="y", alpha=.15)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Parent score − child score\n(common quality-passing models)")
    fig.suptitle("MDM2: one regeneration round, 3 children per arm\nLabels = children meeting all prespecified improvement criteria", fontsize=12)
    fig.savefig(OUT / "analysis/feedback_panel_gain.png", dpi=220)
    fig.savefig(OUT / "analysis/feedback_panel_gain.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
