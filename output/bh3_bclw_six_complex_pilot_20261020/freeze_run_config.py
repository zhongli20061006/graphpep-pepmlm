from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

PILOT = Path(__file__).resolve().parent
ROOT = PILOT.parents[1]
WSL_ROOT = "/mnt/e/Users/graphpep与pepmlm"
PILOT_WSL = f"{WSL_ROOT}/output/bh3_bclw_six_complex_pilot_20261020"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_complex_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence = ""
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, sequence))
            name, sequence = line[1:], ""
        else:
            if name is None:
                raise ValueError("FASTA sequence before header")
            sequence += line
    if name is not None:
        records.append((name, sequence))
    result = []
    for name, sequence in records:
        if sequence.count(":") != 1:
            raise ValueError(f"expected exactly two chains in {name}")
        receptor, peptide = sequence.split(":")
        if len(receptor) != 164 or len(peptide) != 23:
            raise ValueError(f"unexpected chain lengths in {name}")
        result.append((name, receptor, peptide))
    if len(result) != 3:
        raise ValueError("expected exactly three candidate complexes")
    return result


records = read_complex_fasta(PILOT / "complexes.fasta")
peptide_fasta = "".join(f">{name}\n{peptide}\n" for name, _, peptide in records)
(PILOT / "peptide_sequences.fasta").write_text(peptide_fasta, encoding="utf-8")

input_rows = []
for seed in (20261008, 20261009):
    for name, receptor, peptide in records:
        complex_path = (
            f"{PILOT_WSL}/seed_{seed}/relaxed/"
            f"{name}_relaxed_rank_001_alphafold2_multimer_v3_model_1_seed_{seed}.pdb"
        )
        input_rows.append({
            "peptide_id": name,
            "model_id": f"seed{seed}",
            "complex_path": complex_path,
            "peptide_sequence": peptide,
            "receptor_sequence": receptor,
        })
with (PILOT / "graphpep_inputs.csv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(input_rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(input_rows)

foldx_settings = "temperature=298\npH=7.4\nionStrength=0.15\nwater=-IGNORE\nscreen=false\nnoNterm=false\nnoCterm=false\n"
foldx_protocol = {
    "scope": "six relaxed AFM complexes; RepairPDB then AnalyseComplex only; no BuildModel",
    "foldx_executable": "/mnt/e/Users/graphpep与pepmlm/tools/foldx/foldx5_Linux/foldx_20261231",
    "repair_config_template": "command=RepairPDB\npdb=input.pdb\n" + foldx_settings,
    "analyse_config_template": "command=AnalyseComplex\npdb=input_Repair.pdb\nanalyseComplexChains=A,B\n" + foldx_settings,
    "commands": ["foldx -f repair.cfg", "foldx -f analyse.cfg"],
    "score_field": "actual AnalyseComplex A/B Interaction Energy (kcal/mol)",
    "orientation": "lower energy is favorable; report parent energy minus variant energy within seed",
    "failure": "missing/ambiguous Interaction output or command/input error is failed; never substitute BuildModel energy",
}
(PILOT / "foldx_protocol_template.json").write_text(
    json.dumps(foldx_protocol, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)

DRAFT = json.loads((PILOT / "run_config_draft.json").read_text(encoding="utf-8-sig"))
config = DRAFT
config.update({
    "schema": "bclw_six_complex_engineering_pilot_run_config_frozen_v1",
    "created_date": "2026-10-20",
    "status": "FROZEN_PRE_SCORE_NOT_RUN",
    "decision_boundary": (
        "Conditional go only for a six-complex engineering repeatability pilot. "
        "Not an independent affinity benchmark, blind test, or validation of any scorer against Kd."
    ),
    "pre_score_freeze": {
        "frozen_before_any_new_prediction_or_score": True,
        "execution_started": False,
        "execution_authorized_by_this_freeze": False,
        "source_decision": "Dutta 2013 does not add quantitative benchmark pairs; full independent direction benchmark remains on hold.",
        "pilot_value": "Can test local input, structure, mapping, scorer execution, and cross-seed proxy stability only.",
    },
})
config["prediction"]["input_directory"] = PILOT_WSL
config["prediction"]["output_directories"] = [f"{PILOT_WSL}/seed_20261008/structures", f"{PILOT_WSL}/seed_20261009/structures"]
config["prediction"]["num_ensemble"] = 1
config["prediction"]["commands"] = {
    "seed_template": (
        "XLA_PYTHON_CLIENT_PREALLOCATE=false "
        "/home/g/miniconda3/envs/colabfold_exp0/bin/colabfold_batch "
        f"{PILOT_WSL}/complexes.fasta {PILOT_WSL}/seed_<SEED>/structures "
        "--msa-mode single_sequence --model-type alphafold2_multimer_v3 "
        "--num-models 1 --model-order 1 --num-seeds 1 --random-seed <SEED> "
        "--num-recycle 3 --recycle-early-stop-tolerance 0 --num-ensemble 1 "
        "--max-seq 256 --max-extra-seq 512 --num-relax 0 --pair-mode unpaired "
        "--disable-unified-memory --use-fast-kernels --compile-mode fast --sort-queries-by none"
    ),
    "serial_seeds": [20261008, 20261009],
    "remote_msa_or_template_service": False,
    "relax_template": (
        "/home/g/miniconda3/envs/colabfold_exp0/bin/python "
        f"{WSL_ROOT}/pipeline/relax_feedback_panel.py "
        f"{PILOT_WSL}/seed_<SEED>/structures {PILOT_WSL}/seed_<SEED>/relaxed"
    ),
    "relax_settings": {
        "max_iterations": 200,
        "tolerance": 2.39,
        "stiffness": 10.0,
        "max_outer_iterations": 3,
        "use_gpu": False,
    },
}
config["study_scope"]["expected_complex_jobs"] = 6
config["study_scope"]["pairing_rule"] = "Each mutation is compared only with parent under the identical seed; keep both seed-specific differences, do not pool directions."
config["study_scope"]["source_exposure"] = "Public Foight labels were inspected for candidate selection; development-exposed and never blind-test data."
config["candidate_scorers"]["GraphPep"].update({
    "execution_template": (
        f"PATH=/home/g/miniconda3/envs/graphpep_env/bin:$PATH "
        "/home/g/miniconda3/envs/graphpep_env/bin/python "
        f"{WSL_ROOT}/pipeline/step3_graphpep.py --lists {PILOT_WSL}/graphpep_inputs.csv "
        f"--pdb-dir {PILOT_WSL} --scores {PILOT_WSL}/scores/graphpep.csv "
        "--graphpep-dir /home/g/GraphPep_v1.3 "
        f"--patched-resmap {PILOT_WSL}/graphpep_patches/resmap_compat.py --blue-threshold 0.70"
    ),
    "compat_mapper": f"{PILOT_WSL}/graphpep_patches/resmap_compat.py",
    "compat_mapper_sha256": sha(PILOT / "graphpep_patches" / "resmap_compat.py"),
    "pair_metric": "parent raw PRED_SCORE minus variant raw PRED_SCORE, within seed; positive is proxy-only improvement",
})
config["candidate_scorers"]["GraphPep"]["package_code_sha256"]["compat_mapper"] = sha(
    PILOT / "graphpep_patches" / "resmap_compat.py"
)
config["candidate_scorers"]["FoldX"].update({
    "protocol_template": f"{PILOT_WSL}/foldx_protocol_template.json",
    "protocol_template_sha256": sha(PILOT / "foldx_protocol_template.json"),
    "pair_metric": "parent Interaction Energy minus variant Interaction Energy, within seed; positive is proxy-only improvement",
})
config["candidate_scorers"]["ESM2_sequence_prior"].update({
    "implementation_status": "wrapper frozen and helper tests pass; no real checkpoint inference performed",
    "executable": f"{PILOT_WSL}/esm_pll.py",
    "code_sha256": sha(PILOT / "esm_pll.py"),
    "test_file": f"{PILOT_WSL}/test_esm_pll.py",
    "test_sha256": sha(PILOT / "test_esm_pll.py"),
    "command": (
        "PATH=/home/g/miniconda3/envs/graphpep_env/bin:$PATH "
        "/home/g/miniconda3/envs/graphpep_env/bin/python "
        f"{PILOT_WSL}/esm_pll.py {PILOT_WSL}/peptide_sequences.fasta "
        f"{PILOT_WSL}/scores/esm_pll.csv --checkpoint /home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt --device cuda"
    ),
    "output_orientation": "variant sum pseudo-log-likelihood minus parent; higher means sequence-prior preference only",
    "scoring_scope": "one sequence-only score per peptide; not repeated per AFM seed",
})
config["candidate_scorers"]["contacts"].update({
    "no_contact_rule": "zero receptor-peptide residue pairs at <=6.0 A makes the complex a failed row before any scorer",
    "interpretation": "geometry-only; more contacts do not imply stronger affinity",
})
config["proposed_readouts_pre_score"]["quality_gate"]["interface_pae_definition"] = (
    "Maximum PAE over all directed receptor-to-peptide and peptide-to-receptor residue pairs; pass only if <=10 A."
)
config["proposed_readouts_pre_score"]["quality_gate"]["peptide_plddt_definition"] = (
    "Mean C-alpha pLDDT over all peptide-chain residues; pass only if >=70."
)
config["proposed_readouts_pre_score"]["quality_gate"]["iptm_definition"] = "AFM summary ipTM; pass only if >=0.70."
config["proposed_readouts_pre_score"]["quality_gate"]["clash_definition"] = (
    "Minimum heavy-atom interchain distance must be >=1.5 A; otherwise fail the structural gate."
)
config["proposed_readouts_pre_score"]["cross_seed_pose_readout"] = (
    "Report receptor-aligned peptide-backbone RMSD and contact-pair Jaccard descriptively; no pass threshold and no pooled average."
)
config["stop_and_failure_policy"] = [
    "Keep exactly six planned sequence-seed cells in the denominator; do not replace failed outputs, change seeds, or discard failed rows.",
    "Any chain identity/sequence mapping error, quality-gate failure, clash, zero contact, nonfinite or missing score is failed; do not substitute another scorer's value.",
    "A proxy direction is seed-stable only when both seed-specific same-parent comparisons are finite, pass all structural gates, and have the same nonzero sign.",
    "Stop when cumulative accelerator wall time reaches 2 hours or new pilot artifacts exceed 2 GiB; no concurrent GPU jobs.",
    "Do not calculate pooled correlations or claim affinity prediction from two mutation directions.",
    "No PepMLM generation, model training, full-pool prediction, remote MSA request, or P27A minimization is included.",
]
config["budgets"]["two_GiB_interpretation"] = "Hard cap on all new files under this pilot output directory, not GPU memory."
config["budgets"]["execution_policy"] = "Serial local jobs only; no remote services; each new score/prediction counts toward the stated limits."
config["run_authorization"] = {
    "status": "NOT_STARTED",
    "reason": "This turn froze the pre-score configuration only; no explicit instruction to launch the six-complex compute job was given.",
    "pre_run_checks": [
        "Recheck free GPU memory and disk headroom immediately before launch.",
        "Recompute all frozen source/code/checkpoint hashes and verify the six exact PDB paths after prediction.",
        "Launch only if execution is explicitly requested; preserve all failures and raw outputs.",
    ],
}
config["remaining_pre_run_items"] = config["run_authorization"]["pre_run_checks"]
runtime = config["runtime_and_provenance"]
runtime["verified_package_versions"] = {
    "colabfold": "1.6.3",
    "jax": "0.10.2",
    "jaxlib": "0.10.2",
    "openmm": "8.6.1",
    "graphpep_environment_python": "3.11.14",
    "torch": "2.7.1+cu128",
    "fair_esm": "2.0.0",
    "biopython": "1.88",
    "numpy": "1.26.4",
    "pandas": "3.0.5",
}
runtime["command_help_sha256"] = {
    name: sha(PILOT / name)
    for name in (
        "colabfold_batch_help.txt",
        "relax_feedback_panel_help.txt",
        "step3_graphpep_help.txt",
        "esm_pll_help.txt",
    )
}
runtime["run_artifact_sha256"] = {
    "complexes.fasta": sha(PILOT / "complexes.fasta"),
    "peptide_sequences.fasta": sha(PILOT / "peptide_sequences.fasta"),
    "graphpep_inputs.csv": sha(PILOT / "graphpep_inputs.csv"),
    "foldx_protocol_template.json": sha(PILOT / "foldx_protocol_template.json"),
    "esm_pll.py": sha(PILOT / "esm_pll.py"),
    "test_esm_pll.py": sha(PILOT / "test_esm_pll.py"),
    "graphpep_patches/resmap_compat.py": sha(PILOT / "graphpep_patches" / "resmap_compat.py"),
    "pipeline/step3_graphpep.py": sha(ROOT / "pipeline" / "step3_graphpep.py"),
    "pipeline/relax_feedback_panel.py": sha(ROOT / "pipeline" / "relax_feedback_panel.py"),
    "pipeline/make_resmap_compat.py": sha(ROOT / "pipeline" / "make_resmap_compat.py"),
}
runtime["mapper_origin"] = {
    "source": "current workspace pipeline/make_resmap_compat.py BODY",
    "generated_path": f"{PILOT_WSL}/graphpep_patches/resmap_compat.py",
    "generated_sha256": sha(PILOT / "graphpep_patches" / "resmap_compat.py"),
    "legacy_home_copy_sha256_not_used": "114d306b921f9f6b3661b06d85d43eaa849fc7bfb5738225126f0e854acc4fce",
    "note": "The runtime passes this output-local mapper explicitly; the /home/g/pepmlm_pipeline copy is not used as workspace code.",
}
runtime["preflight_scope"] = {
    "cli_help_only": True,
    "real_checkpoint_loaded": False,
    "real_prediction_or_scoring": False,
    "colabfold_flags_checked": [
        "single_sequence", "alphafold2_multimer_v3", "num-models", "model-order",
        "num-seeds", "random-seed", "num-recycle", "recycle-early-stop-tolerance",
        "num-ensemble", "max-seq", "max-extra-seq", "num-relax", "pair-mode",
        "disable-unified-memory", "use-fast-kernels", "compile-mode", "sort-queries-by",
    ],
}
(PILOT / "run_config_frozen_v1.json").write_text(
    json.dumps(config, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8"
)

print(json.dumps({
    "records": len(records),
    "graphpep_rows": len(input_rows),
    "input_fasta_sha256": sha(PILOT / "complexes.fasta"),
    "peptide_fasta_sha256": sha(PILOT / "peptide_sequences.fasta"),
    "graphpep_inputs_sha256": sha(PILOT / "graphpep_inputs.csv"),
    "foldx_template_sha256": sha(PILOT / "foldx_protocol_template.json"),
    "esm_code_sha256": sha(PILOT / "esm_pll.py"),
    "manifest_sha256": sha(PILOT / "run_config_frozen_v1.json"),
}, indent=2))
