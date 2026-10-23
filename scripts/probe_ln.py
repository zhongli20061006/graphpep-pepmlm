#!/usr/bin/env python3
"""List all layer-norm / lm_head / embedding-related keys in the fair-esm checkpoint."""
import torch

CKPT = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"
obj = torch.load(CKPT, map_location="cpu", weights_only=False, mmap=True)
src = obj["model"]
for k in src:
    if ("emb_layer_norm" in k) or ("lm_head" in k) or ("embed" in k) or ("pooler" in k) or ("contact" in k):
        print(f"  {k}  {tuple(src[k].shape)}")
print(f"total {len(src)}")
