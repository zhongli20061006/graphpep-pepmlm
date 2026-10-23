#!/usr/bin/env python3
"""Read fair-esm checkpoint key names + shapes via streaming pickle, low RAM."""
import pickle
import zipfile

CKPT = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"


class Stub:
    def __init__(self, *a, **k):
        self.args = a

    def __repr__(self):
        return f"Stub({self.args[:2]})"


class ProbeUnpickler(pickle.Unpickler):
    def persistent_load(self, pid):
        # Returning a lightweight stub keeps RAM flat; we only want key names here.
        return Stub(*pid)


with zipfile.ZipFile(CKPT) as zf:
    with zf.open("archive/data.pkl") as fh:
        obj = ProbeUnpickler(fh).load()

print("top-level type:", type(obj))
if isinstance(obj, dict):
    print("top-level keys:", list(obj.keys()))
    for k, v in obj.items():
        if isinstance(v, dict):
            ks = list(v.keys())
            print(f"\n=== [{k}] : {len(ks)} keys ===")
            for kk in ks[:12]:
                print("   ", kk)
            print("    ...")
            for kk in ks[-10:]:
                print("   ", kk)
