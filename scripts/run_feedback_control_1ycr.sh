#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

out=output/feedback_control_1ycr_20261006
gp_python=/home/g/miniconda3/envs/graphpep_env/bin/python
cf_python=/home/g/miniconda3/envs/colabfold_exp0/bin/python
afm=/home/g/miniconda3/envs/colabfold_exp0/bin/colabfold_batch
mkdir -p "$out"
exec > >(tee -a "$out/run.log") 2>&1

"$gp_python" pipeline/feedback_control_1ycr.py prepare
XLA_PYTHON_CLIENT_PREALLOCATE=false "$afm" "$out/parent_msa" "$out/parent_structures" \
  --model-type alphafold2_multimer_v3 --num-models 3 --model-order 1,2,3 \
  --num-seeds 1 --random-seed 20261010 --num-recycle 3 --recycle-early-stop-tolerance 0 \
  --max-seq 256 --max-extra-seq 512 --num-relax 0 --disable-unified-memory \
  --use-fast-kernels --compile-mode fast --sort-queries-by none
JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false OPENMM_CPU_THREADS=2 \
  "$cf_python" pipeline/relax_feedback_panel.py "$out/parent_structures" "$out/parent_relaxed"
"$gp_python" pipeline/feedback_control_1ycr.py score-parent
"$gp_python" pipeline/feedback_control_1ycr.py generate
XLA_PYTHON_CLIENT_PREALLOCATE=false "$afm" "$out/offspring_msa" "$out/structures" \
  --model-type alphafold2_multimer_v3 --num-models 3 --model-order 1,2,3 \
  --num-seeds 1 --random-seed 20261010 --num-recycle 3 --recycle-early-stop-tolerance 0 \
  --max-seq 256 --max-extra-seq 512 --num-relax 0 --disable-unified-memory \
  --use-fast-kernels --compile-mode fast --sort-queries-by none
JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false OPENMM_CPU_THREADS=2 \
  "$cf_python" pipeline/relax_feedback_panel.py "$out/structures" "$out/relaxed_structures"
"$gp_python" pipeline/feedback_control_1ycr.py analyze
