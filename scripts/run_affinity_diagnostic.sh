#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
out=output/affinity_diagnostic_mdm2_20261006
gp_python=/home/g/miniconda3/envs/graphpep_env/bin/python
cf_python=/home/g/miniconda3/envs/colabfold_exp0/bin/python
afm=/home/g/miniconda3/envs/colabfold_exp0/bin/colabfold_batch
mkdir -p "$out"
exec > >(tee -a "$out/run.log") 2>&1
"$gp_python" pipeline/affinity_diagnostic.py prepare
for seed in 20261011 20261012; do
  XLA_PYTHON_CLIENT_PREALLOCATE=false "$afm" "$out/afm_msa" "$out/seed_$seed/structures" \
    --model-type alphafold2_multimer_v3 --num-models 3 --model-order 1,2,3 \
    --num-seeds 1 --random-seed "$seed" --num-recycle 3 --recycle-early-stop-tolerance 0 \
    --max-seq 256 --max-extra-seq 512 --num-relax 0 --disable-unified-memory \
    --use-fast-kernels --compile-mode fast --sort-queries-by none
  JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false OPENMM_CPU_THREADS=2 \
    "$cf_python" pipeline/relax_feedback_panel.py "$out/seed_$seed/structures" "$out/seed_$seed/relaxed"
  "$gp_python" pipeline/affinity_diagnostic.py score-afm --seed "$seed"
done
"$gp_python" pipeline/affinity_diagnostic.py score-native
"$gp_python" pipeline/affinity_diagnostic.py analyze
