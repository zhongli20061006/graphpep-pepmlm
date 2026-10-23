#!/usr/bin/env bash
# Install Python deps needed by PepMLM (transformers) and BindEvaluator (sklearn/datasets).
set -u
PIP=/home/g/miniconda3/envs/graphpep_env/bin/pip
echo "=== installing deps $(date -Is) ==="
timeout 2400 $PIP install --disable-pip-version-check -q \
    "transformers>=4.40,<5" \
    scikit-learn \
    datasets \
    seaborn \
    biopython 2>&1 | tail -20
echo "install exit: $?"
echo
echo "=== verify ==="
/home/g/miniconda3/envs/graphpep_env/bin/python - <<'PYEOF'
mods = ["torch", "transformers", "sklearn", "datasets", "pandas", "numpy", "matplotlib",
        "seaborn", "scipy", "esm", "MDAnalysis", "rdkit", "pytorch_lightning"]
for m in mods:
    try:
        mod = __import__(m)
        print(f"  OK   {m:22s} {getattr(mod, '__version__', '?')}")
    except Exception as e:
        print(f"  FAIL {m:22s} {type(e).__name__}: {e}")
PYEOF
echo "=== done $(date -Is) ==="
