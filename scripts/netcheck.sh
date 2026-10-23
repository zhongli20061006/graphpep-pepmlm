#!/usr/bin/env bash
echo "===== A. curl verbose to huggingface ====="
timeout 30 curl -sS -o /dev/null -w "http_code=%{http_code} time=%{time_total}\n" https://huggingface.co 2>&1
echo "curl_exit=$?"
echo
echo "===== B. DNS resolution ====="
timeout 15 getent hosts huggingface.co 2>&1; echo "getent_exit=$?"
timeout 15 getent hosts pypi.org 2>&1; echo "getent_pypi_exit=$?"
echo
echo "===== C. WSL network mode / resolv.conf ====="
cat /etc/resolv.conf 2>&1
echo "--- ip route ---"
ip route 2>&1 | head -5
echo
echo "===== D. ping (may be blocked) ====="
timeout 15 ping -c 2 1.1.1.1 2>&1 | tail -4
echo
echo "===== E. pypi via pip index ====="
timeout 40 /home/g/miniconda3/envs/graphpep_env/bin/pip download --no-deps --dest /tmp/_netcheck idna 2>&1 | tail -5
echo
echo "===== F. proxy env ====="
env | grep -iE "proxy|http_" 2>&1; echo "(proxy end)"
echo
echo "===== G. is git available / git ls-remote ====="
timeout 40 git ls-remote https://huggingface.co/facebook/esm2_t33_650M_UR50D 2>&1 | head -3
echo "git_exit=$?"
