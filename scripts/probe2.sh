#!/usr/bin/env bash
echo "########## T1. BindEvaluator prediction -> binding score mapping ##########"
sed -n '145,182p' /home/g/pepmlm_pipeline/probe/bindevaluator.py
echo
echo "########## T2. bindevaluator_modules/__init__.py exports ##########"
timeout 60 curl -4 -sSL "https://hf-mirror.com/ChatterjeeLab/moPPIt/resolve/main/modules/bindevaluator_modules/__init__.py" 2>&1 | head -40
echo
echo "########## T3. modules.py class signatures ##########"
timeout 60 curl -4 -sSL "https://hf-mirror.com/ChatterjeeLab/moPPIt/resolve/main/modules/bindevaluator_modules/modules.py" -o /tmp/mods.py 2>&1
grep -nE "^class |def __init__|def forward" /tmp/mods.py | head -40
echo
echo "########## T4. models.py ##########"
timeout 60 curl -4 -sSL "https://hf-mirror.com/ChatterjeeLab/moPPIt/resolve/main/modules/bindevaluator_modules/models.py" -o /tmp/mdl.py 2>&1
grep -nE "^class |def __init__|def forward|EsmModel|from_pretrained" /tmp/mdl.py | head -40
echo
echo "########## T5. UniProt P13591 reachable? ##########"
timeout 60 curl -4 -sSL "https://rest.uniprot.org/uniprotkb/P13591.fasta" 2>&1 | head -8
