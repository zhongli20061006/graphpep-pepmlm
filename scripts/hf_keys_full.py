#!/usr/bin/env python3
"""Print the complete HF key mapping target set for the conversion."""
import json
from transformers import EsmConfig, EsmForMaskedLM

with open("/home/g/pepmlm_pipeline/models/PepMLM-650M/config.json") as fh:
    cfg = EsmConfig(**json.load(fh))
sd = EsmForMaskedLM(cfg).state_dict()
keys = list(sd.keys())

print("== top-level groups ==")
groups = {}
for k in keys:
    g = ".".join(k.split(".")[:3])
    groups[g] = groups.get(g, 0) + 1
for g, n in groups.items():
    print(f"   {g}  x{n}")

print("\n== non layer.* keys (everything except esm.encoder.layer.N.*) ==")
for k in keys:
    if not k.startswith("esm.encoder.layer."):
        print(f"   {k}  {tuple(sd[k].shape)}")

print("\n== keys of esm.encoder.layer.0.* (full per-layer template) ==")
for k in keys:
    if k.startswith("esm.encoder.layer.0."):
        print(f"   {k.replace('layer.0.', 'layer.N.')}")

print(f"\ntotal {len(keys)}")
