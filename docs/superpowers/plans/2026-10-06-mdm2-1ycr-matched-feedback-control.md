# MDM2 1YCR Matched Feedback Control Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run a frozen, 1YCR-only comparison of AFM-predicted GraphPep-guided PepMLM regeneration against a matched random-site control, with anchors protected and ten unique offspring per arm.

**Architecture:** Keep the previous multi-parent experiment immutable. Add one focused driver and shell runner that reuse the validated PepMLM, shared receptor MSA, ColabFold AFM-v3, restrained Amber relaxation, and installed GraphPep utilities. Store every new input and result under `output/feedback_control_1ycr_20261006`.

**Tech Stack:** Python, PepMLM-650M in fp16, ColabFold AlphaFold-Multimer v3, CPU Amber minimization, GraphPep at `/home/g/GraphPep_v1.3`, pandas/Biopython.

---

### Task 1: Freeze the controlled experiment protocol

**Files:**
- Create: `pipeline/feedback_control_1ycr.py`
- Modify: `pipeline/feedback_panel.py` (`draw_child` gains optional `top_k=3`; existing callers retain prior behavior)
- Create: `output/feedback_control_1ycr_20261006/protocol.json` (created by the driver before child sequences are scored)

- [ ] **Step 1: Define immutable experiment constants**

The driver must fix these values before generating offspring:

```python
OUT = ROOT / "output/feedback_control_1ycr_20261006"
PARENT_ID = "P53WT_1YCR"
PARENT_SEQUENCE = "ETFSDLWKLLPEN"
PDB = "1YCR"
PROTECTED_POSITIONS = [3, 7, 10]  # F3, W7, L10 in peptide-local numbering
N_CHILDREN_PER_ARM = 10
MUTATIONS_PER_CHILD = 2
TOP_K = 4
TEMPERATURE = 1.0
SEED = 20261010
AFM_MODEL_NUMBERS = [1, 2, 3]
AFM_RECYCLES = 3
```

- [ ] **Step 2: Record site-selection and scoring rules**

Use only parent residues with at least one AFM-parent GraphPep contact observation, then remove positions 3, 7, and 10. The guided arm masks the two remaining positions with the lowest contact-count-weighted GraphPep confidence mean; ties break by ascending peptide position. The random arm draws two distinct positions uniformly from that exact eligible set. Both arms create exactly two substitutions, use PepMLM top-k 4 at temperature 1, and share one global set to prevent duplicate sequences across arms. Top-k 3 would permit only nine joint sequences at the fixed guided pair, so top-k 4 is the smallest equal setting that permits ten unique candidates. Do not filter by PPL.

Reuse the existing structure gate: ipTM >= 0.70, peptide pLDDT >= 70, interface PAE <= 10 Å, minimum inter-chain heavy-atom distance >= 1.5 Å, and at least one interface contact. Compare offspring with the matching AFM parameter model for the parent. A child is an operational GraphPep improvement only if at least two matched models pass, the median parent-minus-child GraphPep score is greater than 0.05, and at least two matched per-model score differences are positive.

- [ ] **Step 3: Save protocol inputs and hashes before PepMLM sampling**

The frozen JSON must include parent sequence and PDB, anchor names and indices, observed eligible sites, guided positions, all random position pairs, candidate counts, all seeds, AFM/relaxation settings, GraphPep source hashes, parent-feedback file hashes, the score gate, and the exact software source hash. Child evaluation must reject a changed protocol or source input.

### Task 2: Prepare parent feedback and generate the two candidate arms

**Files:**
- Modify: `pipeline/feedback_control_1ycr.py`
- Modify: `pipeline/feedback_panel.py`
- Create: `scripts/run_feedback_control_1ycr.sh`
- Create: `output/feedback_control_1ycr_20261006/candidates.csv`
- Create: `output/feedback_control_1ycr_20261006/offspring_msa/*.a3m`

- [ ] **Step 1: Reuse the exact resolved 1YCR parent and validated utilities**

Point the driver to `output/experiment0/expanded_validation/POSCTRL_MDM2_P53_native.pdb`, the cached receptor-only A3M `output/experiment0/expanded_validation/results/POSCTRL_MDM2_P53.a3m`, the local PepMLM model `/home/g/pepmlm_pipeline/models/PepMLM-650M`, and the installed GraphPep `/home/g/GraphPep_v1.3`. Do not download or modify GraphPep.

- [ ] **Step 2: Generate parent AFM structures and freeze feedback before offspring evaluation**

Run AFM-v3 models 1, 2, and 3 on the parent using seed `20261010`, 3 recycles, no templates, max-seq 256, and max-extra-seq 512. Apply the same restrained Amber settings used in the previous experiment to parent and offspring structures. Require at least two parent models to pass the fixed structure gate before creating candidates. Build the predicted parent feedback by pooling contact confidences weighted by contact count across passing parent models.

- [ ] **Step 3: Generate ten guided and ten random children**

Use a separate deterministic seed for every PepMLM draw, exclude the parent amino acid from each masked position's top-k choices, sample the joint masked-position choices by their PepMLM probabilities, and retry through unseen joint alternatives until each arm contains ten unique children. Save exact mask positions, mutations, position-source pool, per-draw seed, selected sequence, and PPL as a diagnostic only. The parent and all 20 offspring use the same receptor-only MSA rows; the query peptide row is updated to each sequence.

### Task 3: Predict, relax, score, and compare all paired structures

**Files:**
- Modify: `pipeline/feedback_control_1ycr.py`
- Modify: `scripts/run_feedback_control_1ycr.sh`
- Create: `output/feedback_control_1ycr_20261006/analysis/model_results.csv`
- Create: `output/feedback_control_1ycr_20261006/analysis/sequence_summary.csv`
- Create: `output/feedback_control_1ycr_20261006/analysis/arm_summary.csv`
- Create: `output/feedback_control_1ycr_20261006/experiment_record.csv`
- Create: `output/feedback_control_1ycr_20261006/report_zh.md`

- [ ] **Step 1: Run the same AFM and Amber settings on parent and offspring**

Predict three AFM parameter models for each of the 20 children plus the parent (63 predicted complexes total). Relax all 63 with identical CPU Amber restraints. Preserve unrelaxed structures and confidence JSON files.

- [ ] **Step 2: Score each relaxed model with the installed GraphPep**

Create 60 offspring GraphPep evaluation records with the locally installed GraphPep code and reuse the three parent scores from this exact AFM seed. A complex with no interface contacts is recorded explicitly as `no_contact`. Match each child model to the parent model with the same AFM model number. Keep raw GraphPep outputs, contact feedback, structure paths, confidence values, and hashes.

- [ ] **Step 3: Produce the paired comparison and integrity record**

Report all ten child-level matched gains per arm, structure-pass counts, improvement counts, median and mean matched gains, mutation positions, and sequence-level notes. Include the protocol and source hashes, verify no protected anchor was mutated, verify exactly two mutations per child and 20 unique child sequences, confirm the receptor MSA rows are identical, and state that the result is a small computational comparison rather than a statistical or experimental-affinity claim.

### Task 4: Execute and review the run

**Files:**
- Verify: `output/feedback_control_1ycr_20261006/report_zh.md`
- Verify: `output/feedback_control_1ycr_20261006/analysis/verification.json`

- [ ] **Step 1: Run the dedicated shell workflow**

Run `bash scripts/run_feedback_control_1ycr.sh` from the project root. Expected completed result: 1 parent plus 20 offspring, 63 relaxed AFM structures, 60 offspring GraphPep evaluation records (including explicit `no_contact` rows), ten complete child records per arm, and a final report.

- [ ] **Step 2: Review all predefined gates and report limitations**

Stop with a recorded incomplete status if fewer than two parent models pass, any child lacks all three AFM outputs, GraphPep scores fail, or protocol/source hashes mismatch. Do not silently change thresholds, seeds, arms, anchors, or site pools after viewing offspring scores.
