# MDM2–p53 feedback pilot: one controlled generation

User approved the PDB-parent pilot on 2026-10-05. Execute in this session using the existing local environments. No training, new model downloads, or native GraphPep edits are required.

## Frozen protocol (before offspring generation)

- Parent: resolved 1YCR MDM2 residues 25–109 (85 aa), p53 residues 17–29 (`ETFSDLWKLLPEN`, 13 aa). This truncated construct is the previously tested positive control, not the full deposited constructs.
- Feedback source: GraphPep on the experimental parent, existing `0015_POSCTRL_MDM2_P53_1YCR_native/residue_feedback.csv`. Eligible positions have at least one contact and a finite mean confidence. Missing contacts are not assigned zero confidence.
- Guided arm: mask the two eligible positions with lowest mean contact confidence, ties broken by position. No literature-based anchor protection is added. Random arm: uniformly sample two positions from the identical eligible set for each replicate. The common eligibility filter means the comparison isolates confidence-based position choice within contacting residues.
- Generate 3 children per arm, single forward pass per proposal, frozen PepMLM, fp16 inference, canonical amino acids only, top-k=3, temperature=1. Exclude the parent amino acid before top-k so every proposal has exactly two substitutions. This is an intentional constrained local-regeneration extension, not verbatim original de novo PepMLM sampling.
- Sample joint two-position alternatives according to the product of masked-position probabilities. Exclude previously sampled identical sequences within each arm to give three distinct proposals; do not filter by PPL or structure score. Seeds: position selection 20261005; categorical proposal draws 20261006–20261008, matched between arms. Save all alternatives, probabilities, substitutions, seeds, and PPL. No cross-arm duplicate exclusion.
- Use the same 85 aa receptor sequence in PepMLM and AFM. Build all seven AFM inputs from the cached receptor-only unpaired MSA; discard peptide homolog and paired rows, replace query peptide, add one peptide-only query. Preserve receptor insertions. Hash the common MSA and verify no parent peptide information persists outside its own query.
- AFM v3 models 1,2,3, seed 20261005, 3 recycles (early-stop tolerance 0), MSA 256/512, no templates, no relaxation, unified memory disabled. Re-predict the parent with the identical settings. Seven sequences × three models = 21 predictions.
- GraphPep scores all models, including poor-quality structures, with quality labels. Lower raw score is better. Score is not affinity. Compare the predicted parent with predicted offspring; native score is a separate diagnostic.
- Predefined exploratory structure gate per model: ipTM >=0.70, peptide mean pLDDT >=70, median bidirectional PAE between contacting receptor residues (6 Å) and peptide <=10 Å, at least one interchain contact, minimum heavy-atom distance >=1.5 Å. These are pilot heuristics, not validated universal thresholds. Sequence gate: at least 2/3 models pass.
- Report all-model median score and spread, quality-pass fraction, quality-gated median, parent-relative changes matched by model, contact-set Jaccard, native contact recall (position-based), peptide CA RMSD after receptor alignment. A failed structure is never ranked as an improved binder merely by score.
- A preliminary improved child must pass the sequence gate, improve median GraphPep score by >0.05 against the predicted parent, and improve in at least 2/3 matched models. Also report strict all-three-model improvement. Compare both arms descriptively; n=3 is not a significance claim.
- PDB feedback initialization uses experimental information; this is local optimization around a known binder. It does not establish prediction-only feedback reliability, de novo generation success, NCAM1 efficacy, binding affinity improvement, or held-out generalization. Native structural similarity is a diagnostic, not proof that alternatives are invalid. Training-set overlap has not been excluded.

## Implementation and verification

- [x] Add `pipeline/feedback_pilot.py`: feedback selection, constrained sampling, input validation, shared-MSA construction, reproducibility metadata.
- [x] Add focused tests for missing-feedback exclusion, insertion-aware A3M handling and peptide-homolog removal, mutation boundaries, and MSA query mismatch rejection. Verify failures before implementation and pass afterward.
- [x] Generate 6 offspring locally; verify all unmasked residues and receptor sequence remain identical and all children have exactly two changes.
- [x] Run all 21 AFM predictions and retain complete logs/configuration.
- [x] Add `pipeline/feedback_pilot_analysis.py`: coordinate/sequence/confidence validation, GraphPep scoring, per-model and per-sequence summaries, fixed-gate comparisons, report.
- [x] Verify chain sequences, model counts, confidence lengths, recomputed GraphPep score formula, source hashes, and native GraphPep source hashes.
- [x] Deliver CSV/JSON records and Chinese report under `output/feedback_pilot_mdm2_20261005` with an explicit positive, negative, or inconclusive result.
