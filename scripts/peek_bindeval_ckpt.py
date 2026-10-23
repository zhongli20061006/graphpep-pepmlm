#!/usr/bin/env python3
"""Inspect tensor names/shapes inside the partially downloaded BindEvaluator checkpoint.

A torch>=2.6 checkpoint is a ZIP; its central directory sits at the END of the file. We
already have a non-contiguous prefix in parts/, plus we can range-read the tail from the
mirror, so we can enumerate every tensor name and dtype WITHOUT the full 2.8 GB file.

Purpose: confirm our rebuilt PeptideModel will `load_state_dict(..., strict=True)` cleanly,
        and read off the trained hyper-parameters, before the download even finishes.
"""
from __future__ import annotations

import collections
import io
import os
import struct
import subprocess
import sys
import urllib.request
import zipfile

D = "/home/g/pepmlm_pipeline/models/moPPIt/classifier_ckpt"
PARTS = os.path.join(D, "ranges")
TOTAL = 2832618887
URL = ("https://hf-mirror.com/ChatterjeeLab/moPPIt/resolve/main/"
       "classifier_ckpt/finetuned_BindEvaluator.ckpt")
CACHE = "/tmp/bindeval_tail.zip"


def fetch_range(start: int, end: int, out: str) -> None:
    cmd = ["curl", "-4", "-sSL", "--max-time", "900", "-r", f"{start}-{end}", "-o", out, URL]
    subprocess.run(cmd, check=True)


def main() -> int:
    have = 0
    if os.path.isdir(PARTS):
        have = sum(os.path.getsize(os.path.join(PARTS, f))
                   for f in os.listdir(PARTS) if f.startswith("part_"))
    print(f"downloaded so far: {have/1e6:.1f} MB of {TOTAL/1e6:.1f} MB")

    # ZIP End-Of-Central-Directory: scan the final 128 KB
    tail_len = min(131072, TOTAL)
    print(f"range-reading the final {tail_len} bytes ...")
    if not os.path.isfile(CACHE) or os.path.getsize(CACHE) != tail_len:
        fetch_range(TOTAL - tail_len, TOTAL - 1, CACHE)
    tail = open(CACHE, "rb").read()
    print(f"got tail: {len(tail)} bytes")

    # locate EOCD (or ZIP64 EOCD locator)
    idx = tail.rfind(b"PK\x05\x06")
    if idx < 0:
        print("no EOCD found in tail; trying a bigger tail")
        return 1
    cd_size, cd_offset = struct.unpack("<II", tail[idx + 12: idx + 20])
    print(f"central directory: offset={cd_offset} size={cd_size}")

    # For ZIP64 the real values live in the ZIP64 EOCD record
    z64 = tail.rfind(b"PK\x06\x06")
    if z64 >= 0:
        cd_size64, cd_offset64 = struct.unpack("<QQ", tail[z64 + 40: z64 + 56])
        print(f"ZIP64 central directory: offset={cd_offset64} size={cd_size64}")
        cd_offset, cd_size = cd_offset64, cd_size64

    print(f"range-reading central directory ({cd_size} bytes) ...")
    cd_path = "/tmp/bindeval_cd.bin"
    fetch_range(cd_offset, cd_offset + cd_size - 1, cd_path)
    cd = open(cd_path, "rb").read()
    print(f"got central directory: {len(cd)} bytes")

    # parse central directory entries
    names = []
    p = 0
    while p + 46 <= len(cd) and cd[p:p + 4] == b"PK\x01\x02":
        comp = struct.unpack("<H", cd[p + 10:p + 12])[0]
        csize = struct.unpack("<I", cd[p + 20:p + 24])[0]
        usize = struct.unpack("<I", cd[p + 24:p + 28])[0]
        nlen = struct.unpack("<H", cd[p + 28:p + 30])[0]
        elen = struct.unpack("<H", cd[p + 30:p + 32])[0]
        clen = struct.unpack("<H", cd[p + 32:p + 34])[0]
        name = cd[p + 46:p + 46 + nlen].decode("utf-8", "replace")
        names.append((name, csize, usize))
        p += 46 + nlen + elen + clen
    print(f"\nzip members: {len(names)}")
    for n, c, u in names[:12]:
        print(f"  {n}  comp={c}  uncomp={u}")

    # the state_dict pickle is archive/data.pkl -- stream it with a stub persistent_load
    pkl = [n for n, _, _ in names if n.endswith("data.pkl")]
    if not pkl:
        print("no data.pkl member found")
        return 1
    pkl_name = pkl[0]
    entry = [e for e in names if e[0] == pkl_name][0]
    hdr_off = None
    # local file header offset is at cd entry +42
    p = 0
    while p + 46 <= len(cd) and cd[p:p + 4] == b"PK\x01\x02":
        nlen = struct.unpack("<H", cd[p + 28:p + 30])[0]
        elen = struct.unpack("<H", cd[p + 30:p + 32])[0]
        clen = struct.unpack("<H", cd[p + 32:p + 34])[0]
        name = cd[p + 46:p + 46 + nlen].decode("utf-8", "replace")
        if name == pkl_name:
            hdr_off = struct.unpack("<I", cd[p + 42:p + 46])[0]
            csize = struct.unpack("<I", cd[p + 20:p + 24])[0]
            usize = struct.unpack("<I", cd[p + 24:p + 28])[0]
            break
        p += 46 + nlen + elen + clen

    print(f"\n{pkl_name}: local header offset={hdr_off} comp={csize} uncomp={usize}")
    if hdr_off is None or usize > 200_000_000:
        print("pkl too large to range-read comfortably; using torch.load on what we have instead")
        return 0

    print(f"range-reading the pickle ({usize} bytes) ...")
    pkl_path = "/tmp/bindeval_data.pkl.bin"
    fetch_range(hdr_off, hdr_off + csize - 1, pkl_path)
    raw = open(pkl_path, "rb").read()

    # strip the local file header to get at the compressed payload
    sig, ver, flags, method, mtime, mdate, crc, csz, usz, nlen, elen = struct.unpack(
        "<IHHHHHIIIHH", raw[:30])
    payload = raw[30 + nlen + elen:]
    print(f"  method={method} nlen={nlen} elen={elen} payload={len(payload)}")

    if method == 8:
        import zlib
        data = zlib.decompress(payload, -15)
    else:
        data = payload
    print(f"  decompressed pickle: {len(data)} bytes")

    import pickle

    class Stub:
        def __init__(self, *a):
            self.args = a

    class PU(pickle.Unpickler):
        def persistent_load(self, pid):
            return Stub(*pid)

    try:
        obj = PU(io.BytesIO(data)).load()
    except Exception as exc:  # noqa: BLE001
        print(f"  pickle load stopped early: {type(exc).__name__}: {exc}")
        return 0

    print("\n=== recovered object ===")
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, dict):
                print(f"[{k}] dict with {len(v)} entries")
                if k in ("state_dict", "model"):
                    ks = list(v.keys())
                    print("  first 18 keys:")
                    for kk in ks[:18]:
                        shp = tuple(v[kk].shape) if hasattr(v[kk], "shape") else "?"
                        print(f"    {kk}  {shp}")
                    print("  ...")
                    for kk in ks[-8:]:
                        shp = tuple(v[kk].shape) if hasattr(v[kk], "shape") else "?"
                        print(f"    {kk}  {shp}")
            else:
                print(f"[{k}] = {str(v)[:120]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
