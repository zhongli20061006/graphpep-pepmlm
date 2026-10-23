# Experiment 0 cloud pilot implementation plan

> **For agentic workers:** Execute this approved pilot inline, task by task. No cloud charges or model training are authorized.

**Goal:** Prepare two reproducible NCAM1–peptide cloud predictions and a verified local path for importing and scoring their real structures.

**Architecture:** Preserve the installed GraphPep v1.3. Produce an embedded-input Colab notebook from the official batch notebook. Import each predicted complex separately, retaining its own receptor and peptide coordinates. Save confidence and contact feedback with explicit model IDs; raw GraphPep scores rank ascending.

**Tech Stack:** Existing WSL Python 3.11 environment, standard library, pandas, Biopython, installed GraphPep, official ColabFold AFM v3 cloud runtime.

## Task 1: Regression checks

- [x] Add `tests/test_experiment0.py` using unittest and real synthetic PDB/CSV fixtures.
- [x] Run regression checks in the existing WSL environment and observe failures for shared-receptor inputs, missing pilot preparation, and legacy report behavior.
- [x] Verify pair coordinates, preserved sidechains, deterministic candidate selection, two-chain FASTA, full input provenance, score direction, and missing-contact handling.

## Task 2: Prepare the pilot

- [x] Create `pipeline/experiment0.py` with `prepare` and `import` subcommands.
- [x] Deduplicate by sequence; select minimum-PPL sequence and the most distant sequence in the middle PPL quintile, breaking ties by composition difference and proximity to pool median.
- [x] Verify NCAM1 input equals current UniProt P13591 standard sequence before slicing positions 20–718. Record source checksum, receptor checksum, range, and assumption of extracellular binding.
- [x] Produce pilot CSV, two-chain FASTA, receptor FASTA, config, notebook, Chinese instructions, and pending results template.
- [x] Adapt the pinned official notebook. Embed real input, use five AFM v3 models, three recycles, one seed, no relaxation, and save timing/GPU/parameters in a downloadable archive.

## Task 3: Correct local structure scoring

- [x] Replace shared-receptor/mixed-sequence batching with one input pair per manifest row.
- [x] Support cloud paths/model IDs; validate sequences and preserve heavy atoms and coordinates.
- [x] Precheck six-Angstrom contacts and leave empty graphs unscored.
- [x] Preserve raw feedback, join explicit IDs, and restore original residue numbering.
- [x] Correct report sorting/comparisons to lower-is-better and exclude missing feedback from all-blue claims.

## Task 4: Import and verify

- [x] Implement/test import of AFM PDBs and matching confidence JSON, including swapped chain order and canonical residue mapping.
- [x] Run six regression checks, notebook Python/bash syntax validation, and real-data preparation.
- [x] Verify author decoy 1: revised score -3.5795, author reference -3.5795; 58 contact pairs. Keep control separate from candidates.
- [x] Synchronize changed scripts with backups and verify installed GraphPep source hashes remain unchanged.
- [x] Deliver the notebook and instructions. Cloud prediction remains awaiting user-started Colab runtime because this session has no logged-in browser control tool.

## Current state

Two candidates: PEP00815 / TKKRQAARGK (PPL 6.4646) and PEP00356 / TIQKKLLGIP (PPL 19.3150), Hamming distance 9.
Canonical P13591 sequence verified; receptor 699 aa, each complex 709 aa. No candidate cloud structure or candidate GraphPep score exists yet.
