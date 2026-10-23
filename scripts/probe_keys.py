#!/usr/bin/env python3
"""Print fair-esm checkpoint key structure WITHOUT loading all tensors into RAM."""
import zipfile
import sys

CKPT = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"

print(f"opening zip: {CKPT}")
with zipfile.ZipFile(CKPT) as zf:
    names = zf.namelist()
    print(f"zip entries: {len(names)}")
    # peek the top-level data.pkl to list keys
    pkl = [n for n in names if n.endswith("data.pkl")]
    print("pkl entries:", pkl)
    import pickle
    with zf.open(pkl[0]) as fh:
        obj = pickle.load(fh)
    print(f"top-level object type: {type(obj)}")
    if isinstance(obj, dict):
        print("top-level keys:", list(obj.keys())[:20])
        for k, v in obj.items():
            if isinstance(v, dict):
                ks = list(v.keys())
                print(f"\n[{k}] {len(ks)} keys; first 15:")
                for kk in ks[:15]:
                    print("   ", kk)
                print("   ... last 8:")
                for kk in ks[-8:]:
                    print("   ", kk)
