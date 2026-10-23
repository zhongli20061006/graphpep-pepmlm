#!/usr/bin/env python3
"""Ground truth: what does fair-esm itself predict for the same masked position?

If fair-esm also says 'L', my HF conversion is fine and my test expectation was wrong.
If fair-esm says 'Q', the conversion is broken.
"""
import torch
import esm

CKPT = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"

print("loading fair-esm alphabet/model ...", flush=True)
model, alphabet = esm.pretrained.load_model_and_alphabet_local(CKPT)
model.eval()
batch_converter = alphabet.get_batch_converter()

seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYG"
masked = seq[:10] + "<mask>" + seq[11:]
print("masked seq:", masked, " true char@11:", seq[10])

data = [("p", masked)]
_, _, tokens = batch_converter(data)
print("tokens shape:", tuple(tokens.shape))
mask_idx = (tokens[0] == alphabet.mask_idx).nonzero().flatten().tolist()
print("mask token position(s):", mask_idx, "alphabet.mask_idx:", alphabet.mask_idx)

with torch.no_grad():
    out = model(tokens)
logits = out["logits"][0, mask_idx[0]]

tok_to_idx = alphabet.tok_to_idx
aa_order = [(t, tok_to_idx[t]) for t in alphabet.all_toks[:33]]
print("\ntop10 predictions (fair-esm raw):")
top = torch.topk(logits, 10)
for v, i in zip(top.values, top.indices):
    print(f"   {alphabet.get_tok(int(i))}  {float(v):+.3f}")

true_tok = seq[10]
true_i = tok_to_idx[true_tok]
print(f"\ntrue token '{true_tok}' logit={float(logits[true_i]):+.3f} "
      f"rank={(logits > logits[true_i]).sum().item()}")
