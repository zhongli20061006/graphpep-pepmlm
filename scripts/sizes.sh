#!/usr/bin/env bash
export HF_ENDPOINT=https://hf-mirror.com
echo "########## MODEL FILE SIZES (via hf-mirror) ##########"
for repo in ChatterjeeLab/PepMLM-650M ChatterjeeLab/moPPIt ChatterjeeLab/moPPIt-v2; do
  echo "=== $repo ==="
  timeout 45 curl -4 -sS "https://hf-mirror.com/api/models/$repo?blobs=true" 2>&1 \
    | /home/g/miniconda3/envs/graphpep_env/bin/python -c "
import sys, json
try:
    d = json.load(sys.stdin)
except Exception as e:
    print('  parse-fail', e); sys.exit()
print('  downloads:', d.get('downloads'), 'lib:', d.get('library_name'), 'pipeline:', d.get('pipeline_tag'))
tot = 0
for s in d.get('siblings', []):
    sz = s.get('size')
    if sz:
        tot += sz
        print(f\"  {s['rfilename']:<55} {sz/1e6:10.1f} MB\")
    else:
        print(f\"  {s['rfilename']:<55} {'?':>10}\")
print(f'  TOTAL ~ {tot/1e9:.2f} GB')
"
  echo
done

echo "########## ESMFOLD weights on mirror ##########"
for u in \
  "https://hf-mirror.com/facebook/esmfold_v1/resolve/main/config.json" \
  "https://hf-mirror.com/facebook/esmfold_v1/resolve/main/pytorch_model.bin" \
  "https://hf-mirror.com/facebook/esm2_t33_650M_UR50D/resolve/main/pytorch_model.bin" ; do
  timeout 40 curl -4 -sSI "$u" 2>&1 | grep -iE "^HTTP|content-length" | head -2 | sed "s|^|  $(basename $(dirname $u))/$(basename $u): |"
done
