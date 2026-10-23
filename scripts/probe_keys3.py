#!/usr/bin/env python3
"""Inspect fair-esm checkpoint keys using mmap (low RAM)."""
import gc
import sys

import torch

CKPT = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"
print("loading with mmap=True ...", flush=True)
obj = torch.load(CKPT, map_location="cpu", weights_only=False, mmap=True)
print("top-level type:", type(obj))
if isinstance(obj, dict):
    print("top-level keys:", list(obj.keys()))
    for k, v in obj.items():
        if isinstance(v, dict):
            ks = list(v.keys())
            print(f"\n=== [{k}] : {len(ks)} keys ===")
            for kk in ks[:10]:
                shp = tuple(v[kk].shape) if hasattr(v[kk], "shape") else "?"
                print(f"    {kk}  {shp}")
            print("    ...")
            for kk in ks[-6:]:
                shp = tuple(v[kk].shape) if hasattr(v[kk], "shape") else "?"
                print(f"    {kk}  {shp}")
del obj
gc.collect()
print("\nok")
