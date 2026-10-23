#!/usr/bin/env python3
"""Definitive equivalence test: converted HF model vs original fair-esm model.

Feeds identical token id sequences to both and compares the LM logits.
If they match numerically, the conversion is exact and any "wrong" masked-token
prediction is simply what the real ESM-2 650M predicts.
"""
import argparse
import esm
import torch

CKPT = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"
LOCAL = "/home/g/pepmlm_pipeline/models/esm2_t33_650M_UR50D"

# torch>=2.6 defaults weights_only=True; the checkpoint legitimately contains Namespace.
_orig_load = torch.load


def patched_load(*a, **kw):
    kw.setdefault("weights_only", False)
    return _orig_load(*a, **kw)


esm.pretrained.torch.load = patched_load

print("=== loading original fair-esm model ===", flush=True)
fair_model, alphabet = esm.pretrained.load_model_and_alphabet_local(CKPT)
fair_model.eval()

print("=== loading converted HF model ===", flush=True)
from transformers import AutoTokenizer, EsmForMaskedLM

tok = AutoTokenizer.from_pretrained(LOCAL)
hf_model = EsmForMaskedLM.from_pretrained(LOCAL).eval()

seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYG"
masked = seq[:10] + "<mask>" + seq[11:]

# fair-esm representation (batched: (1, L))
bc = alphabet.get_batch_converter()
_, _, ftokens = bc([("p", masked)])
# HF representation
htok = tok(masked, return_tensors="pt")["input_ids"]
print("fair tokens :", ftokens[0].tolist())
print("hf   tokens :", htok[0].tolist())
assert ftokens[0].tolist() == htok[0].tolist(), "tokenisation differs!"

with torch.no_grad():
    f_logits = fair_model(ftokens)["logits"][0]          # (L, 33)
    h_logits = hf_model(htok)["logits"][0]               # (L, 33)

print("\nshapes:", tuple(f_logits.shape), tuple(h_logits.shape))
diff = (f_logits - h_logits).abs()
print(f"max |diff| = {diff.max().item():.6e}")
print(f"mean|diff| = {diff.mean().item():.6e}")

mask_pos = (ftokens[0] == alphabet.mask_idx).nonzero().flatten().tolist()
print("mask positions:", mask_pos)
for p in mask_pos:
    ft = torch.topk(f_logits[p], 5)
    ht = torch.topk(h_logits[p], 5)
    print(f"\npos {p}: true token = '{masked[p]}'")
    print("  fair-esm top5:", [(alphabet.get_tok(int(i)), round(float(v), 3)) for v, i in zip(ft.values, ft.indices)])
    print("  HF       top5:", [(tok.convert_ids_to_tokens([int(i)])[0], round(float(v), 3)) for v, i in zip(ht.values, ht.indices)])

ok = diff.max().item() < 1e-3
print(f"\nEQUIVALENCE: {'PASS' if ok else 'FAIL'} (max diff {diff.max().item():.3e})")
