#!/usr/bin/env bash
# Fetch the moPPIt source files BindEvaluator needs (code only, small).
set -u
M="https://hf-mirror.com/ChatterjeeLab/moPPIt/resolve/main"
D=/home/g/pepmlm_pipeline/models/moPPIt
mkdir -p "$D"
FILES="
bindevaluator.py
moppit.py
README.md
scripts/bindevaluator.sh
scripts/moppit.sh
modules/__init__.py
modules/general_module.py
modules/bindevaluator_modules/__init__.py
modules/bindevaluator_modules/dataloaders.py
modules/bindevaluator_modules/layers.py
modules/bindevaluator_modules/models.py
modules/bindevaluator_modules/modules.py
modules/bindevaluator_modules/score_domain.py
models/__init__.py
models/classifier.py
models/peptide_classifiers.py
models/peptide_models.py
models/peptiverse_classifiers.py
utils/__init__.py
utils/dataloader.py
utils/dataset.py
utils/flow_utils.py
utils/parsing.py
classifier_code/binding_affinity_unpooled.py
classifier_code/__init__.py
"
for f in $FILES; do
  out="$D/$f"
  mkdir -p "$(dirname "$out")"
  code=$(curl -4 -sSL -o "$out" -w "%{http_code}" --retry 3 "$M/$f" 2>/dev/null)
  sz=$(stat -c '%s' "$out" 2>/dev/null || echo 0)
  printf "  %-58s code=%s size=%s\n" "$f" "$code" "$sz"
done
echo
echo "=== tree ==="
find "$D" -type f | sort
