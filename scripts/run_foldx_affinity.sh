#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

foldx_output=output/foldx_affinity_mdm2_20261006
foldx_python=/home/g/miniconda3/envs/graphpep_env/bin/python
foldx_binary="${1:-${FOLDX_BINARY:-}}"
mkdir -p "$foldx_output"
exec > >(tee -a "$foldx_output/run.log") 2>&1

"$foldx_python" pipeline/foldx_affinity.py prepare
"$foldx_python" pipeline/verify_foldx_inputs.py
if [[ -z "$foldx_binary" ]]; then
  printf '%s\n' 'FoldX 尚未提供：输入已准备；实际能量评分未运行。' >&2
  printf '%s\n' '用法：bash scripts/run_foldx_affinity.sh /absolute/path/to/authorized/foldx' >&2
  exit 3
fi

"$foldx_python" pipeline/foldx_affinity.py run --foldx "$foldx_binary"
"$foldx_python" pipeline/foldx_affinity.py analyze
"$foldx_python" pipeline/foldx_affinity.py verify
"$foldx_python" pipeline/verify_foldx_scores.py
