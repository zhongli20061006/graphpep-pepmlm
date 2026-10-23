#!/usr/bin/env python3
"""Step 4 - compare the two filtering strategies and write output/stage1_report.md.

The primary, statistically clean comparison available without the structure stage is on
the 1000-candidate pool itself: how do the PPL ranking and the BindEvaluator ranking agree,
and where do their top sets land relative to each other?

When output/graphpep_scores.csv exists (structure stage complete) the report additionally
carries the GraphPep Top-10 comparison table, box/violin figure and conclusion.
Until then it is emitted as an explicitly-marked interim report.

Definitions used:
  * ppl_score        PepMLM pseudo-perplexity (authors' compute_pseudo_perplexity); lower is better
  * bind_eval_score  BindEvaluator site-probabilities averaged over residues with p >= 0.5
                     (the threshold in the authors' own main()); higher is better, 0.0 if none
  * binding_score    GraphPep PRED_SCORE = -log(1 + pred_fnat); lower is more confident
  * red_residue_ratio fraction of contact pairs with confidence < BLUE_THRESHOLD;
                     "all-blue" means the ratio is exactly 0
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/g/pepmlm_pipeline"
OUT = os.path.join(ROOT, "output")
BLUE_THRESHOLD = 0.70


def read_id_list(path: str) -> list[str]:
    if not os.path.isfile(path):
        return []
    ids = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("peptide_id"):
                continue
            ids.append(line.split("\t")[0])
    return ids


def group_stats(df: pd.DataFrame, col: str) -> dict:
    s = pd.to_numeric(df[col], errors="coerce").dropna()
    return {
        "n": int(len(s)),
        "mean": float(s.mean()) if len(s) else float("nan"),
        "median": float(s.median()) if len(s) else float("nan"),
        "std": float(s.std(ddof=1)) if len(s) > 1 else 0.0,
        "min": float(s.min()) if len(s) else float("nan"),
        "max": float(s.max()) if len(s) else float("nan"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", default=os.path.join(OUT, "generated_peptides.csv"))
    ap.add_argument("--top-ppl", default=os.path.join(OUT, "top20_ppl.txt"))
    ap.add_argument("--top-bindeval", default=os.path.join(OUT, "top20_bindeval.txt"))
    ap.add_argument("--scores", default=os.path.join(OUT, "graphpep_scores.csv"))
    ap.add_argument("--report", default=os.path.join(OUT, "stage1_report.md"))
    ap.add_argument("--figure", default=os.path.join(OUT, "boxplot_comparison.png"))
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--force-interim", action="store_true")
    args = ap.parse_args()

    has_graphpep = os.path.isfile(args.scores) and not args.force_interim
    pep = pd.read_csv(args.peptides) if os.path.isfile(args.peptides) else None
    ppl_ids = read_id_list(args.top_ppl)
    be_ids = read_id_list(args.top_bindeval)
    overlap = sorted(set(ppl_ids) & set(be_ids))

    L: list[str] = []
    A = L.append
    A("# Stage 1 report — PepMLM generation → BindEvaluator filtering → GraphPep rescoring")
    A("")
    A("**Target** NCAM1 / CD56 (UniProt `P13591`), 858 aa &nbsp;·&nbsp; "
      "**Generator** PepMLM-650M &nbsp;·&nbsp; **Filters** PPL vs BindEvaluator (moPPIt) "
      "&nbsp;·&nbsp; **Downstream scorer** GraphPep v1.3")
    A("")
    A(f"Status: {'**complete**' if has_graphpep else '**interim — structure stage pending**'}")
    A("")

    # ------------------------------------------------------------------ 1
    A("## 1. Experiment summary")
    A("")
    if pep is None:
        A("- candidate table not found; nothing to summarise.")
        A("")
    else:
        uniq = pep["peptide_sequence"].nunique()
        lens = sorted(pep["peptide_sequence"].str.len().unique().tolist())
        A(f"- candidates generated: **{len(pep)}** (unique sequences **{uniq}**)")
        A(f"- peptide lengths present: {lens} (requested length 10, top-k = 3)")
        st = group_stats(pep, "ppl_score")
        A(f"- `ppl_score`: min {st['min']:.3f} · median {st['median']:.3f} · "
          f"max {st['max']:.3f} (lower = more target-consistent)")
        if "bind_eval_score" in pep.columns:
            sb = group_stats(pep, "bind_eval_score")
            A(f"- `bind_eval_score`: min {sb['min']:.4f} · median {sb['median']:.4f} · "
              f"max {sb['max']:.4f}")
            n0 = int((pep["bind_eval_score"] == 0).sum())
            A(f"- candidates whose predicted binding-site set was **empty** (score 0.0): "
              f"**{n0}/{len(pep)}** ({100*n0/len(pep):.1f}%)")
        A(f"- Top-20 lists: PPL {len(ppl_ids)} ids, BindEvaluator {len(be_ids)} ids; "
          f"**overlap {len(overlap)}** peptides")

        # candidate-pool quality (from pool_stats.py, if present)
        ps_path = os.path.join(OUT, "pool_stats.json")
        if os.path.isfile(ps_path):
            import json as _json
            ps = _json.load(open(ps_path))
            A("")
            A("### Candidate-pool quality (diversity and composition)")
            A("")
            A(f"- duplicate fraction: **{ps['duplicate_fraction']*100:.1f}%** "
              f"({ps['n_candidates']} candidates, {ps['n_unique']} unique)")
            A(f"- mean per-position amino-acid entropy (0–1, 1 = uniform over 20 aa): "
              f"**{ps['mean_positional_entropy']:.3f}**")
            A(f"- per-position consensus (majority residue and its share): "
              f"{', '.join(ps['position_consensus'])}")
            top_pool = sorted(ps["composition_pct"].items(), key=lambda kv: -kv[1])[:5]
            top100 = sorted(ps["composition_pct_top100_ppl"].items(), key=lambda kv: -kv[1])[:5]
            A(f"- pool composition (top 5): " +
              ", ".join(f"{aa} {p:.1f}%" for aa, p in top_pool))
            A(f"- top-100 by PPL composition (top 5): " +
              ", ".join(f"{aa} {p:.1f}%" for aa, p in top100))
            A("")
            A("The pool is **strongly compositionally biased** and the mean positional "
              "entropy is low: with `top-k = 3` at every position, sampling concentrates on "
              "the few residues the model already prefers for this target. Selecting the "
              "lowest-PPL candidates makes this *worse*, not better "
              "(the PPL-top-100 is more K-rich than the pool). A PPL-only filter therefore "
              "trades diversity for target-consistency, which is exactly the limitation a "
              "local mask-resampling loop is meant to address.")
            A("")
    A("")

    # ------------------------------------------------------------------ 2
    A("## 2. Metric definitions")
    A("")
    A("| metric | definition | direction |")
    A("|---|---|---|")
    A("| `ppl_score` | PepMLM pseudo-perplexity of the binder region — mask one binder "
      "position at a time, take the true token's cross-entropy, average, `exp()`. "
      "Reproduces the authors' `compute_pseudo_perplexity`. | lower better |")
    A("| `bind_eval_score` | BindEvaluator's per-residue sigmoid binding-site probabilities "
      "on the target, averaged over residues with `p >= 0.5` — the threshold hard-coded in "
      "the authors' `main()`. `0.0` when no residue passes. | higher better |")
    A(f"| `binding_score` | GraphPep `PRED_SCORE = -log(1 + pred_fnat)` "
      "(`bin/core.py::_score`); a pose/contact-quality score. | lower better |")
    A(f"| `red_residue_ratio` | fraction of a peptide's contact pairs whose GraphPep "
      f"confidence is `< {BLUE_THRESHOLD:.2f}`; \"all-blue\" = ratio exactly 0. | lower better |")
    A("")
    A("> BindEvaluator is a **sequence-only** model. It does not output a single peptide-level "
      "binding probability, so `bind_eval_score` is an explicit aggregation choice "
      "(documented above) rather than a model-native score. "
      "`output/bindeval_diagnostics.csv` carries `be_n_sites`, `be_max_prob` and "
      "`be_mean_prob` so the choice can be re-derived or replaced without re-running the model.")
    A("")

    # ------------------------------------------------------------------ 3
    A("## 3. Do the two filters select the same candidates?")
    A("")
    if pep is None or "bind_eval_score" not in pep.columns:
        A("bind_eval_score not available yet.")
        A("")
    else:
        sub = pep.dropna(subset=["ppl_score", "bind_eval_score"])
        if sub["bind_eval_score"].nunique() <= 1:
            A("> **`bind_eval_score` is constant across all candidates**, so it cannot "
              "discriminate and the comparison is not informative. See the note in section 6.")
            A("")
        r_pear = sub["ppl_score"].corr(sub["bind_eval_score"], method="pearson")
        r_spear = sub["ppl_score"].corr(sub["bind_eval_score"], method="spearman")
        A(f"- Pearson r(PPL, bind_eval_score) = **{r_pear:+.4f}**")
        A(f"- Spearman ρ(PPL, bind_eval_score) = **{r_spear:+.4f}**")
        A(f"- overlap of the two Top-20 lists = **{len(overlap)}/20**")
        A("")
        A("| ranking | top-5 candidates (id, sequence, score) |")
        A("|---|---|")
        top_ppl = sub.nsmallest(5, "ppl_score")
        top_be = sub.nlargest(5, "bind_eval_score")
        A("| by PPL (lowest) | " + "; ".join(
            f"`{r.peptide_id}` {r.peptide_sequence} ({r.ppl_score:.2f})"
            for r in top_ppl.itertuples()) + " |")
        A("| by bind_eval_score (highest) | " + "; ".join(
            f"`{r.peptide_id}` {r.peptide_sequence} ({r.bind_eval_score:.3f})"
            for r in top_be.itertuples()) + " |")
        A("")
        A(f"Interpretation: a near-zero correlation with low Top-20 overlap means the two "
          f"filters are selecting largely **different** candidates — the choice of filter "
          f"materially changes which peptides advance to the structure stage.")
        A("")

        # top-10 comparison on each filter's own scale
        A(f"### Top-{args.top_n} score summaries (on the filter's own scale)")
        A("")
        n = args.top_n
        s_p = group_stats(sub.nsmallest(n, "ppl_score"), "ppl_score")
        s_b = group_stats(sub.nlargest(n, "bind_eval_score"), "bind_eval_score")
        A("| group | n | mean | median | std | min | max |")
        A("|---|---|---|---|---|---|---|")
        A(f"| PPL-filtered Top-{n} (`ppl_score`) | {s_p['n']} | {s_p['mean']:.4f} | "
          f"{s_p['median']:.4f} | {s_p['std']:.4f} | {s_p['min']:.4f} | {s_p['max']:.4f} |")
        A(f"| BindEvaluator Top-{n} (`bind_eval_score`) | {s_b['n']} | {s_b['mean']:.4f} | "
          f"{s_b['median']:.4f} | {s_b['std']:.4f} | {s_b['min']:.4f} | {s_b['max']:.4f} |")
        A("")

        # figure: PPL vs bind_eval_score scatter + Top-20 set diagram
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
        ax = axes[0]
        ax.scatter(sub["ppl_score"], sub["bind_eval_score"], s=14, alpha=0.45,
                   color="#4C72B0", label=f"all {len(sub)} candidates")
        ax.scatter(sub.nsmallest(20, "ppl_score")["ppl_score"],
                   sub.nsmallest(20, "ppl_score")["bind_eval_score"],
                   s=42, color="#C44E52", label="Top-20 by PPL", zorder=3)
        ax.set_xlabel("PepMLM pseudo-perplexity (lower = better)")
        ax.set_ylabel("bind_eval_score (higher = better)")
        ax.set_title(f"Filter disagreement: r={r_pear:+.3f}, ρ={r_spear:+.3f}")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

        ax = axes[1]
        only_ppl = len(set(ppl_ids) - set(be_ids))
        only_be = len(set(be_ids) - set(ppl_ids))
        span = len(set(ppl_ids) | set(be_ids))
        ax.barh(["PPL only", "shared", "BindEvaluator only"], [only_ppl, len(overlap), only_be],
                color=["#4C72B0", "#8C8C8C", "#DD8452"])
        for i, v in enumerate([only_ppl, len(overlap), only_be]):
            ax.text(v + 0.15, i, str(v), va="center", fontsize=10)
        ax.set_title(f"Top-20 membership (union = {span})")
        ax.set_xlabel("peptides")
        ax.grid(axis="x", alpha=0.3)
        plt.rcParams["axes.unicode_minus"] = False
        fig.suptitle("PPL vs BindEvaluator filtering on the shared candidate pool")
        fig.tight_layout()
        fig.savefig(args.figure, dpi=200)
        print(f"figure written: {args.figure}")
        A(f"![Filter comparison]({os.path.basename(args.figure)})")
        A("")

    # ------------------------------------------------------------------ 4
    A("## 4. GraphPep structure stage")
    A("")
    if not has_graphpep:
        A("> **Not yet run.** `output/graphpep_scores.csv` does not exist, so the "
          "GraphPep binding-score comparison, the box/violin figure and the "
          "which-filter-is-better conclusion are intentionally **absent** from this report "
          "rather than guessed.")
        A("")
        A("Why it is deferred, and what it needs:")
        A("")
        A("1. **ESMFold v1 weights are 8.44 GB** and were not present on this machine; at the "
          "observed mirror throughput that is a multi-hour download.")
        A("2. More importantly, **ESMFold cannot predict peptide-protein complexes.** It is a "
          "single-chain predictor with no interface-confidence output, so the peptide's "
          "placement against the receptor is not a binding pose.")
        A("3. The task's comparison needs *peptide-dependent* geometry. Any deterministic "
          "placement (fixed offset, surface tangent, hotspot anchoring) produces "
          "**near-identical scores for every peptide**, which would make the two groups "
          "statistically indistinguishable and the conclusion vacuous.")
        A("4. NCAM1 offers no template interface to fall back on: all 7 PDB entries for "
          "P13591 are isolated Ig/FN3 domains (`5AEA` Ig1 1.9 Å, `3MTR` Ig5-FN1 1.8 Å, "
          "`2VKW`/`2VKX` FN3 1-2, `2E3V`/`2HAZ`/`5LKN` FN3 1) — none contains a bound peptide.")
        A("")
        A("The active route is cloud AlphaFold-Multimer v3: prepare the pilot with "
          "`pipeline/experiment0.py prepare`, import real complexes with "
          "`pipeline/experiment0.py import`, then score each receptor/peptide pair "
          "separately with `pipeline/step3_graphpep.py`. Pilot results live in "
          "`output/experiment0/` and do not complete the Top-20 filter comparison.")
        A("")
    else:
        sc = pd.read_csv(args.scores)
        sc["binding_score"] = pd.to_numeric(sc["binding_score"], errors="coerce")
        sc = sc[np.isfinite(sc["binding_score"])]
        if "n_contacts" in sc.columns:
            sc = sc[pd.to_numeric(sc["n_contacts"], errors="coerce") > 0]
        sc["source"] = sc["source"].astype(str).str.strip().str.lower()
        g_ppl = sc[sc["source"] == "ppl"]
        g_be = sc[sc["source"].isin(["bindeval", "bind_eval", "be"])]
        n = args.top_n
        t_ppl = g_ppl.nsmallest(n, "binding_score")
        t_be = g_be.nsmallest(n, "binding_score")

        def gp_stats(d: pd.DataFrame) -> dict:
            s = pd.to_numeric(d["binding_score"], errors="coerce").dropna()
            rr = pd.to_numeric(d["red_residue_ratio"], errors="coerce").dropna()
            return {
                "n": int(len(s)),
                "mean": float(s.mean()), "median": float(s.median()),
                "std": float(s.std(ddof=1)) if len(s) > 1 else 0.0,
                "min": float(s.min()), "max": float(s.max()),
                "all_blue_n": int((rr == 0).sum()),
                "all_blue_valid_n": len(rr),
                "all_blue_frac": float((rr == 0).mean()) if len(rr) else float("nan"),
            }

        s_p, s_b = gp_stats(t_ppl), gp_stats(t_be)
        A(f"### Top-{n} comparison on the GraphPep binding score")
        A("")
        A(f"| group | n | mean | median | std | min | max | all-blue (ratio 0) |")
        A("|---|---|---|---|---|---|---|---|")
        for name, s in (("PPL-filtered", s_p), ("BindEvaluator-filtered", s_b)):
            A(f"| {name} | {s['n']} | {s['mean']:.4f} | {s['median']:.4f} | {s['std']:.4f} | "
              f"{s['min']:.4f} | {s['max']:.4f} | {s['all_blue_n']}/{s['all_blue_valid_n']} "
              f"({s['all_blue_frac']*100:.1f}%) |")
        A("")
        try:
            from scipy.stats import mannwhitneyu
            if len(t_ppl) >= 2 and len(t_be) >= 2:
                u, p = mannwhitneyu(t_ppl["binding_score"], t_be["binding_score"],
                                    alternative="two-sided")
                A(f"Mann-Whitney U = {u:.1f}, p = {p:.4f} — **exploratory only**; the two "
                  f"groups come from one candidate pool and are not independent.")
                A("")
        except Exception as exc:  # noqa: BLE001
            print(f"mann-whitney skipped: {exc}")

        # box/violin figure requested by the task
        fig2, axes2 = plt.subplots(1, 2, figsize=(11, 4.6))
        data = [t_ppl["binding_score"].astype(float).values,
                t_be["binding_score"].astype(float).values]
        labels = ["PPL-filtered", "BindEvaluator-filtered"]
        colors = ["#4C72B0", "#DD8452"]
        ax = axes2[0]
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False)
        ax.set_xticks([1, 2], labels)
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.65)
        rng = np.random.default_rng(0)
        for i, d in enumerate(data, start=1):
            ax.scatter(np.full(len(d), i) + rng.normal(0, 0.045, len(d)), d, s=26,
                       color="black", alpha=0.7, zorder=3)
        ax.set_ylabel("GraphPep pose/contact score (lower = better)")
        ax.set_title(f"Top-{n} binding score")
        ax.grid(axis="y", alpha=0.3)
        ax = axes2[1]
        nonconstant = [(i, d) for i, d in enumerate(data, start=1) if len(d) > 1 and np.ptp(d) > 0]
        parts = ax.violinplot([d for _, d in nonconstant], positions=[i for i, _ in nonconstant],
                              showmeans=True, showextrema=False) if nonconstant else {"bodies": []}
        for body, c in zip(parts["bodies"], [colors[i-1] for i, _ in nonconstant]):
            body.set_facecolor(c)
            body.set_alpha(0.6)
        for i, d in enumerate(data, start=1):
            if len(d):
                ax.scatter(np.full(len(d), i), d, color=colors[i-1], s=22)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(labels)
        ax.set_ylabel("GraphPep binding score")
        ax.set_title(f"Top-{n} distribution")
        ax.grid(axis="y", alpha=0.3)
        fig2.suptitle("GraphPep binding score: PPL filtering vs BindEvaluator filtering")
        fig2.tight_layout()
        gp_fig = args.figure
        fig2.savefig(gp_fig, dpi=200)
        print(f"graphpep figure written: {gp_fig}")
        A(f"![Binding score comparison]({os.path.basename(gp_fig)})")
        A("")
        diff = s_b["mean"] - s_p["mean"]
        better = "BindEvaluator" if diff < 0 else "PPL"
        A("### Conclusion")
        A("")
        A(f"- Mean Top-{n} GraphPep binding score: PPL-filtered {s_p['mean']:.4f} vs "
          f"BindEvaluator-filtered {s_b['mean']:.4f} (Δ = {diff:+.4f}).")
        if np.isfinite(diff) and diff != 0:
            A(f"- **{better} filtering retained candidates with the lower mean GraphPep score** "
              f"under this protocol.")
        else:
            A("- The available scores do not establish a difference between filters.")
        A(f"- All-blue fraction: PPL {s_p['all_blue_frac']*100:.1f}% vs "
          f"BindEvaluator {s_b['all_blue_frac']*100:.1f}%.")
        A("")

    # ------------------------------------------------------------------ 5
    A("## 5. Limitations")
    A("")
    A("- PepMLM generation here follows the authors' **published Space implementation** "
      "(single forward pass, top-k sampling at mask positions). It is not an iterative "
      "masked-refinement loop.")
    A("- PPL was computed with a batched fast path, verified against the authors' literal "
      "per-peptide loop (relative deviation reported in `logs/`).")
    A("- `bind_eval_score` is an aggregation choice, not a model-native score (see section 2).")
    A("- GraphPep scores are only meaningful for complexes with real inter-chain contact; its "
      "contact threshold is 6 Å (`bin/config.py`).")
    A("- The two Top lists are small (n=10–20) and derived from the same candidate pool, so "
      "they are **not independent samples**; any p-value would be exploratory.")
    A("")

    os.makedirs(os.path.dirname(args.report), exist_ok=True)
    with open(args.report, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"report written: {args.report}  ({len(L)} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
