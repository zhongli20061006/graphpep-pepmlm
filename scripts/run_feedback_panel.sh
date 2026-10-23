#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
panel_out=output/feedback_panel_mdm2_20261005
panel_gp_python=/home/g/miniconda3/envs/graphpep_env/bin/python
panel_cf_python=/home/g/miniconda3/envs/colabfold_exp0/bin/python
panel_afm=/home/g/miniconda3/envs/colabfold_exp0/bin/colabfold_batch
exec > >(tee "$panel_out/panel_console.log") 2>&1
"$panel_gp_python" pipeline/feedback_panel.py score-parents
"$panel_gp_python" pipeline/feedback_panel.py generate
XLA_PYTHON_CLIENT_PREALLOCATE=false "$panel_afm" "$panel_out/offspring_msa" "$panel_out/structures" \
  --model-type alphafold2_multimer_v3 --num-models 3 --model-order 1,2,3 \
  --num-seeds 1 --random-seed 20261009 --num-recycle 3 --recycle-early-stop-tolerance 0 \
  --max-seq 256 --max-extra-seq 512 --num-relax 0 --disable-unified-memory \
  --use-fast-kernels --compile-mode fast --sort-queries-by none
JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false OPENMM_CPU_THREADS=2 \
  "$panel_cf_python" pipeline/relax_feedback_panel.py "$panel_out/structures" "$panel_out/relaxed_structures"
"$panel_gp_python" pipeline/feedback_panel.py analyze
