# Fixed versus flexible peptide-backbone diagnostic implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans. The user explicitly approved execution. Continue through implementation, checks, 48-structure pilot and full 240-structure run without another permission request. This directory is not a Git repository; no worktree, commit or merge step applies.

**Goal:** Test whether allowing local peptide-backbone movement changes the previously failed FoldX P27A direction while retaining the L22A weaker-binding control.

**Architecture:** Leave all prior experiments and scripts unchanged. Load 60 frozen WT/mutant pairs (two mutations, six backbones, five repeats). Prepare one deterministic hydrogenated coordinate set for each of 120 input structures; reuse it across fixed and flexible conditions. Minimize with existing OpenMM CPU and score with the licensed FoldX. Expand from preselected repeat 1 to all five repeats on technical validity, irrespective of affinity-score direction. Independently audit raw files and aggregates.

**Tech Stack:** Existing WSL `/home/g/miniconda3/envs/colabfold_exp0/bin/python`, OpenMM 8.6.1 CPU, NumPy/Biopython, Python standard-library state/logging, official FoldX 5.1. No training, new AFM inference, dependency downloads, MD or experimental Kd measurement.

## Frozen protocol and ownership

- Scientific settings and 164 workspace plus six external hashes are already frozen in `output/backbone_diagnostic_20261007/protocol.json` before new relaxation or FoldX scores.
- Create `pipeline/backbone_diagnostic.py`: CLI and immutable inputs, common hydration, phase selection, isolated workers, command logs, tracked results, resume, source/resource verification.
- Create `pipeline/backbone_relax.py`: OpenMM preparation/minimization and exact constraints; no FoldX or scientific outcome selection.
- Create `pipeline/backbone_geometry.py`: validated coordinates, displacement, torsions, descriptive alpha-like geometry, O--N distances and receptor-aligned comparisons.
- Create `pipeline/backbone_analysis.py`: raw-derived structure/pair/backbone/seed summaries, corresponding flexible-minus-fixed contrasts, actual energy terms and a Chinese report; no importing scores as truth.
- Create `tests/test_backbone_diagnostic.py`: meaningful scientific/input/state/geometry guards, fixtures in temporary directories.
- Root owns `output/backbone_diagnostic_20261007/independent_audit.py` and final supplemental interpretation/figures. Do not edit prior frozen files, protocol or environment.

## Tasks

- [x] Inspect installed OpenMM/platform/XML interfaces; freeze protocol before producing new scientific outputs.
- [x] Write tests first, run RED, implement the four modules, then obtain GREEN. Check frozen versus mobile atoms, units (0.1 kcal/mol/A² = 41.84 kJ/mol/nm²), paired gain sign, median of paired contrasts, missing results, malformed raw tables, changed provenance, and incomplete versus complete phase states. Root fresh verification: 15 tests passed in 23.243 s.
- [x] Independently review scientific/spec compliance, then code quality; fix substantive findings. Freeze implementation hashes after reviewers approve and edits stop. Spec PASS; quality APPROVE after completed-but-untracked hydration recovery and residue-count units were corrected.
- [x] Verify an isolated real preparation/minimization probe. Confirm hydrogen preparation preserves heavy-atom coordinates and chemical identity; confirm actual fixed atoms remain unchanged, flexible peptide moves, finite energies and correct serialization. Probe results remain excluded from formal statistics. Actual max-2000 probe: fixed 223 reporter calls / 74.32 s / accepted force RMS, flexible 224 calls / 73.80 s / not accepted under the stricter free-vector RMS; frozen displacement 0 in both and flexible B max 1.0625 A. No FoldX probe scores.
- [x] Publish baseline diagnostics of all 18 independently predicted repaired WT/P27A/L22A AFM structures and all 60 prior pairs' energy-term differences. Retain receptor-aligned peptide RMSD, per-residue phi/psi, local B8--11 geometry and O_i--N_i+4 distances. Independent checks passed: 18 structures, 216 torsion rows, 144 distance-proxy rows and 1,740 prior term contrasts (count terms use count units). These are geometric proxies, not measured hydrogen-bond energies or DSSP assignments.
- [x] Run the 48 preselected repeat-1 structure jobs in two isolated CPU workers. Preserve raw outputs. Record technical pilot gate independently of desired score direction. Completed 48/48 valid structures and 24/24 corresponding pairs; technical PASS, 6/48 meet the stricter free-vector RMS acceptance measure. Full phase started automatically at 09:54 China time; scientific score direction did not select expansion.
- [x] On technical validity, run remaining repeats through all 240 structure jobs / 120 condition-specific pairs. Pilot files remain the corresponding preselected subset; never double-count them. Completed 240/240 valid structures, 120/120 paired rows and 60/60 paired contrasts; no failed jobs.
- [x] Analyze and verify all results; independently reparse actual PDBs/tables and recompute paired gains, contrasts, five-run then three-model summaries and signs. Check all frozen input, implementation, resource and published-artifact hashes. Independent audit passed: 164 frozen workspace sources, six external sources, 480 raw tables and 8,160 numeric terms checked; maximum frozen-atom displacement 0 Å.
- [x] Publish exact counts, convergence/constraint outcomes, uncertainty, figures and interpretation; preserve old continuation gate false. Results: 35/240 met the strict free-atom force-RMS threshold; P27A direction was not recovered in either full seed; L22A controls were correct in both. The descriptive local alpha-like count increased under flexibility in both P27A mutant and WT reference structures, with median mutation-specific change 0; O7···N11 remained above the 3.5 Å distance-proxy threshold in all 30 mutant structures. Reports and figure are in `output/backbone_diagnostic_20261007/analysis/`.

## Scientific constraints

- Input identities: A is the existing full 85-aa MDM2 sequence, B is 12-mer p53 WT `ETFSDLWKLLPE`, P27A `ETFSDLWKLLAE` or L22A `ETFSDAWKLLPE`; preserve numbering 1..85 / 1..12 and both OXT termini.
- Fixed condition sets A/B N, CA, C, O, OXT particle masses to zero. Flexible condition sets only A backbone masses to zero, restores normal B masses and applies a weak B-backbone coordinate tether. All other degrees of freedom receive identical handling. Do not impose helical torsions.
- Add hydrogens once per input role with deterministic preparation seed; use the same coordinates in both arms. Input heavy-atom coordinates must remain unchanged during hydration.
- Freeze ff14SB/OBC2, dielectric values, kappa, CPU deterministic one-thread settings, 10 kJ/mol/nm minimization tolerance and 2,000-iteration cap as recorded in protocol. Measure finite energies, actual iterations and force RMS; report nonconvergence honestly rather than treating the iteration cap as convergence.
- Strip only H/D before FoldX; no further RepairPDB. Both conditions use identical AnalyseComplex A/B settings and actual corresponding WT reference. Require exactly one Interaction and Summary table with agreeing finite energies. Do not use OpenMM potential energy or FoldX folding stability as binding free energy.
- Full aggregation requires five corresponding gains per backbone and three complete backbone medians per seed. Compute flexible-minus-fixed difference within each same repeat before aggregation, not by independently subtracting marginal medians.
- Failure/provenance states remain explicit, auditable and non-promotable; a failed job must not create a valid paired gain. Do not selectively omit a negative or failed result.
- Technical pilot gate: all 48 jobs identity/terminal/finite/raw-table valid, frozen atoms unchanged, at least one flexible peptide displacement above 0.01 Å; report convergence. No gate on biological score sign.
- Fixed-versus-flexible differences support only local geometric sensitivity. The weak tether, implicit solvent, fixed receptor, lack of unbound/thermal sampling and forcefield error limit interpretation. Seeds/models/repeats are not biological replicates. Historical 1YCR is a different 13-mer and cannot silently replace the 12-mer assay construct.

## Run interfaces

```bash
cd /mnt/e/Users/graphpep与pepmlm
/home/g/miniconda3/envs/colabfold_exp0/bin/python -m unittest discover -s tests -p test_backbone_diagnostic.py -v
/home/g/miniconda3/envs/colabfold_exp0/bin/python pipeline/backbone_diagnostic.py prepare
/home/g/miniconda3/envs/colabfold_exp0/bin/python pipeline/backbone_diagnostic.py baseline
/home/g/miniconda3/envs/colabfold_exp0/bin/python pipeline/backbone_diagnostic.py run --phase pilot --workers 2
/home/g/miniconda3/envs/colabfold_exp0/bin/python pipeline/backbone_diagnostic.py run --phase full --workers 2
/home/g/miniconda3/envs/colabfold_exp0/bin/python pipeline/backbone_diagnostic.py analyze
/home/g/miniconda3/envs/colabfold_exp0/bin/python pipeline/backbone_diagnostic.py verify
```

Expected final counts: 120 prepared/common-hydrogenated inputs, 240 minimized structures and successful AnalyseComplex outputs, 120 paired rows, 60 corresponding flexible-minus-fixed contrasts, 24 condition/backbone summaries, eight condition/seed summaries.

## Sources

- [OpenMM User Guide](https://docs.openmm.org/latest/userguide/application/02_running_sims.html): implicit solvent, forcefields, local minimization and platforms.
- [OpenMM minimizer tests](https://github.com/openmm/openmm/blob/master/tests/TestLocalEnergyMinimizer.h): frozen/mass-zero particle behavior, verified again in the actual installed engine.
- [Li et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/): reference assay and proposed conformational contributions; this remains a previously inspected development panel.
