#!/usr/bin/env bash
# Read-only reconnaissance for the pepmlm_pipeline baseline experiment.
echo "===== 1. GraphPep models/ ====="
ls -la /home/g/GraphPep_v1.3/models/ 2>&1
echo
echo "===== 2. GraphPep ckpt/ ====="
ls -la /home/g/GraphPep_v1.3/ckpt/ 2>&1
echo
echo "===== 3. Search for PepMLM / moPPIt / ESMFold ====="
find /home/g -maxdepth 5 \( -iname "*pepmlm*" -o -iname "*moppit*" -o -iname "*esmfold*" \) 2>/dev/null | head -30
echo "(search end)"
echo
echo "===== 4. HuggingFace cache ====="
ls -la /home/g/.cache/huggingface 2>&1 | head -20
du -sh /home/g/.cache/huggingface 2>/dev/null
echo
echo "===== 5. Torch cache ====="
ls -la /home/g/.cache/torch 2>&1 | head -20
echo
echo "===== 6. Internet check ====="
timeout 25 curl -s -o /dev/null -w "huggingface.co -> %{http_code}\n" https://huggingface.co 2>&1
timeout 25 curl -s -o /dev/null -w "hf esm2 api -> %{http_code}\n" https://huggingface.co/api/models/facebook/esm2_t33_650M_UR50D 2>&1
echo
echo "===== 7. Disk free ====="
df -h /home/g /tmp | tail -3
echo
echo "===== 8. tmux sessions ====="
tmux list-sessions 2>&1
echo
echo "===== 9. Running heavy procs ====="
ps -ef | grep -E "predict[.]py|run_case|train[.]py" | grep -v grep | head
echo "(proc end)"
echo
echo "===== 10. graphpep_env versions ====="
/home/g/miniconda3/envs/graphpep_env/bin/python -c "import torch, sys; print('python', sys.version.split()[0]); print('torch', torch.__version__); print('cuda_avail', torch.cuda.is_available()); print('dev', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'n/a')" 2>&1 | tail -6
echo
echo "===== 11. installed pkgs of interest ====="
/home/g/miniconda3/envs/graphpep_env/bin/pip list 2>/dev/null | grep -iE "torch|esm|transformers|lightning|geometric|rdkit|mdanalysis|peft|accelerate|bitsandbytes|numpy|pandas|matplotlib|scipy|biopython|openbabel" 2>&1
echo
echo "===== 12. GraphPep.sh key settings ====="
grep -nE "GraphPep_root|esm2_path|conda|CUDA|python" /home/g/GraphPep_v1.3/GraphPep.sh 2>&1 | head -30
