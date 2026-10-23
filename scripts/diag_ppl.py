#!/usr/bin/env python3
"""Diagnose the fast-path vs reference pseudo-PPL disagreement."""
import torch
import torch.nn.functional as F
from transformers import AutoModelForMaskedLM, AutoTokenizer

MODEL_DIR = "/home/g/pepmlm_pipeline/models/PepMLM-650M"
device = "cuda" if torch.cuda.is_available() else "cpu"

tok = AutoTokenizer.from_pretrained(MODEL_DIR)
model = AutoModelForMaskedLM.from_pretrained(MODEL_DIR).to(device).eval()

seq = "".join(l.strip() for l in open("/home/g/pepmlm_pipeline/data/targets/NCAM1.fasta")
              if not l.startswith(">")).upper()
pep = "GIKQLAIIIK"
L = len(pep)

# ---- 1. tokenization equivalence ----
a = tok(seq + pep)["input_ids"]
b = tok(seq)["input_ids"] + tok(pep)["input_ids"][1:-1]
print("tokenize(concat) len      =", len(a))
print("tokenize(target)+pep body =", len(b))
print("identical                 =", a == b)
if a != b:
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            print(f"  first diff at {i}: {x} ({tok.convert_ids_to_tokens([x])}) "
                  f"vs {y} ({tok.convert_ids_to_tokens([y])})")
            break
print("tail of a:", a[-14:], tok.convert_ids_to_tokens(a[-14:]))
print("tail of b:", b[-14:], tok.convert_ids_to_tokens(b[-14:]))

# ---- 2. loss extraction equivalence ----
tensor_input = tok.encode(seq + pep, return_tensors="pt").to(device)
Lfull = tensor_input.shape[1]
print("\nfull len:", Lfull, " n_pos:", L)

ref_losses = []
with torch.no_grad():
    for i in range(-L - 1, -1):
        mi = tensor_input.clone()
        mi[0, i] = tok.mask_token_id
        labels = torch.full(tensor_input.shape, -100).to(device)
        labels[0, i] = tensor_input[0, i]
        ref_losses.append(model(mi, labels=labels).loss.item())
print("reference per-position losses:", [round(x, 5) for x in ref_losses])

fast_losses = []
with torch.no_grad():
    for j in range(L):
        row = tensor_input.clone()
        row[0, Lfull - L + j] = tok.mask_token_id
        logits = model(row).logits
        logp = F.log_softmax(logits[0, Lfull - L + j, :].float(), dim=-1)
        tgt = tensor_input[0, Lfull - L + j]
        fast_losses.append(-logp[tgt].item())
print("fast      per-position losses:", [round(x, 5) for x in fast_losses])

print("\nper-position diffs:", [round(f - r, 6) for f, r in zip(fast_losses, ref_losses)])
import numpy as np
print(f"\nREF  exp(mean) = {np.exp(np.mean(ref_losses)):.6f}")
print(f"FAST exp(mean) = {np.exp(np.mean(fast_losses)):.6f}")

# ---- 3. does the attention mask matter? ----
print("\n--- check: is <eos> position handled the same? ---")
print("tensor_input tail:", tensor_input[0, -13:].tolist())
print("the reference loop masks positions -L-1..-1 =",
      list(range(Lfull - L - 1, Lfull)))
