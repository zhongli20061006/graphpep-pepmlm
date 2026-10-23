# FoldX matched WT/mutant diagnostic implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. The user has authorized execution in this conversation; proceed through the checks without another approval request.

**Goal:** Determine whether local FoldX mutation/repacking estimates are more consistent than independently predicted AFM structures for four experimentally characterized MDM2 peptide mutations.

**Architecture:** Preserve the completed 60-model affinity screen unchanged. In a new output directory, mutate each of six existing repaired WT structures per family with FoldX BuildModel, five runs per mutation, then AnalyseComplex each generated mutant and its own corresponding WT reference. Audit raw pair identities, sequences, backbone displacement, energy signs and aggregation before interpreting the diagnostic.

**Tech Stack:** Existing WSL graphpep_env Python 3.11, standard-library CSV/JSON/subprocess, licensed FoldX 5.1. No new AFM inference, model training or downloads.

## Fixed scientific protocol

- Output: `output/foldx_matched_diagnostic_20261006`.
- Parents: full 12-mer p53 `ETFSDLWKLLPE` and PMI `TSFAEYWNLLSP`, receptor A chain is the existing 85-aa MDM2 25–109 sequence. Peptide is chain B, numbered 1–12; free N/C termini retained.
- Mutations, in execution order: p53 P27A = `PB11A;` (stronger), p53 L22A = `LB6A;` (weaker), PMI N8A = `NB8A;` (stronger), PMI W7A = `WB7A;` (weaker).
- Sources: both existing seeds 20261011/20261012, AFM parameter models 1/2/3; only the twelve WT `input_Repair.pdb` files, no selecting a preferred source model after scoring.
- Each of 24 jobs: BuildModel `numberOfRuns=5`, `moveNeighbours=true`, `out-pdb=true`. Use PdbList to identify the five specific mutant/WT pairs, verifying identities rather than assuming a shared WT.
- Fixed settings: 298 K, pH 7.4, ionStrength 0.15 M, water=-IGNORE, noNterm=false, noCterm=false. Keep official defaults for other parameters and log actual commands/outputs. Five repeats assess local optimization reproducibility, not biological replication or unbound peptide conformational sampling.
- Metric: AnalyseComplex(A,B) Interaction Energy. Paired gain = E(WT-reference) − E(mutant), kcal/mol. Never treat BuildModel Dif Total Energy as binding ΔΔG.
- Record full AnalyseComplex numeric terms, including both intraclashes columns. Record backbone displacement from the starting structure and within each pair; preserve terminal atoms and require complete sequences.
- Aggregation: median five gains within each source backbone; then median three source-backbone medians per seed, retaining all individual gains/ranges/sign counts. Missing/failed pairs remain explicitly missing; no selective drop or pass declaration from partial results.
- Diagnostic label per mutation/seed: correct only if all three complete backbone medians agree with the experimental direction; mixed or zero is inconclusive. This is a new diagnostic, not a changed affinity-screen gate or external validation.
- Compare to original independently predicted AFM matched gains, which remain unchanged. Improvements under a common starting structure suggest sensitivity to geometry/repacking; they do not prove AFM causation or reliable affinity prediction. A fixed-backbone approach can miss helix/conformational effects of P27A/N8A.

## Tasks and checks

- [x] Freeze `protocol.json` with mutation identities, settings, all input/reference/binary/resource SHA-256 values before the isolated real-format probe. Probe outputs stay under `runtime_preflight` and are excluded from formal results.
- [x] Add `tests/test_foldx_matched_diagnostic.py` before implementation; run it to demonstrate missing-feature failure. Fixtures stay in TemporaryDirectory and are explicitly synthetic. Test ambiguous/nonexistent PdbList pairs, incorrect source/chain/residue numbering, wrong mutation and missing terminal atoms, gain sign, median-of-paired-differences, failed/partial jobs, malformed/nonfinite energy tables and altered frozen hashes.
- [x] Create `pipeline/foldx_matched_diagnostic.py`: prepare 24 immutable jobs; parse PdbList; run BuildModel and paired AnalyseComplex; preserve command logs and artifact hashes; resume only verified complete jobs; analyze raw data into pair/backbone/seed summaries and energy-term CSVs plus Chinese report.
- [x] Probe the actual official BuildModel/PdbList output with p53 P27A and five runs; confirm complete WT/mutant identities and pair count before interpreting energies. No scientific parameter tuning from the probe.
- [x] Review specification, then code quality; resolve substantive issues before the full run. Freeze implementation hashes once tests and reviews pass.
- [x] Run all 24 jobs, at most two isolated CPU jobs concurrently. Publish progress and failures; preserve interrupted attempts rather than deleting them.
- [x] Independently parse every actual output, recalculate five-run and three-backbone medians and signs, compare published CSVs, verify source/tool/artifact hashes and reject incomplete results.
- [x] Update the new README and report with actual counts and uncertainty; leave the original screen protocol, decisions, structures and scripts unchanged.

Completed 2026-10-06: 23/23 tests; 24 successful jobs, 120 complete pairs, 240 successful AnalyseComplex calls; independent audit passed (480 raw tables, 8,160 numeric terms; backbone displacement 0 Å). P27A remained incorrect in both seeds; N8A and both weaker controls had correct signs. Old continuation gate remains false; no new experimental Kd or holdout validation.

## Commands

```bash
cd /mnt/e/Users/graphpep与pepmlm
/home/g/miniconda3/envs/graphpep_env/bin/python -m unittest discover -s tests -p test_foldx_matched_diagnostic.py -v
/home/g/miniconda3/envs/graphpep_env/bin/python pipeline/foldx_matched_diagnostic.py prepare
/home/g/miniconda3/envs/graphpep_env/bin/python pipeline/foldx_matched_diagnostic.py run --foldx /mnt/e/Users/graphpep与pepmlm/tools/foldx/foldx5_Linux/foldx_20261231
/home/g/miniconda3/envs/graphpep_env/bin/python pipeline/foldx_matched_diagnostic.py analyze
/home/g/miniconda3/envs/graphpep_env/bin/python pipeline/foldx_matched_diagnostic.py verify
```

Expected complete counts: 24 BuildModel jobs, 120 matched pairs, 240 AnalyseComplex outputs. Scientific results can be positive, negative or inconclusive; technical correctness never substitutes for affinity validation.

## Sources

- [FoldX BuildModel](https://foldxsuite.crg.eu/command/BuildModel): matching WT controls and five-repeat convergence checks.
- [Individual-list mutation syntax](https://foldxsuite.crg.eu/parameter/mutant-file).
- [FoldX AnalyseComplex](https://foldxsuite.crg.eu/command/AnalyseComplex): Interaction Energy and component tables.
- [Li et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/): same-assay reference affinities and conformational mechanisms; already-used development panel, not holdout.
