#!/usr/bin/env python3
"""Build EsmForMaskedLM from config (random init) and dump its exact state_dict keys.

This gives the ground-truth key names for the fair-esm -> transformers conversion,
without needing any network access.
"""
import json
import torch
from transformers import EsmConfig, EsmForMaskedLM

CFG = "/home/g/pepmlm_pipeline/models/PepMLM-650M/config.json"
with open(CFG) as fh:
    cfgd = json.load(fh)

cfg = EsmConfig(**cfgd)
model = EsmForMaskedLM(cfg)
sd = model.state_dict()

print(f"HF EsmForMaskedLM: {len(sd)} keys")
keys = list(sd.keys())
for k in keys[:22]:
    print(f"   {k}  {tuple(sd[k].shape)}")
print("   ...")
for k in keys[-14:]:
    print(f"   {k}  {tuple(sd[k].shape)}")
