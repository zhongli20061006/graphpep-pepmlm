#!/usr/bin/env python3
"""Profile PepMLM forward-pass throughput: fp32 vs fp16, and batch scaling.

Goal: find the cheapest safe configuration for the 1000-peptide run on an 8 GB GPU.
"""
import time
import torch
import torch.nn.functional as F
from transformers import AutoModelForMaskedLM, AutoTokenizer

MD = "/home/g/pepmlm_pipeline/models/PepMLM-650M"
device = "cuda"
tok = AutoTokenizer.from_pretrained(MD)
seq = "".join(l.strip() for l in open("/home/g/pepmlm_pipeline/data/targets/NCAM1.fasta")
              if not l.startswith(">")).upper()
ids1 = tok(seq + "A" * 10)["input_ids"]
print("seq tokens:", len(ids1))

for dtype_name, dtype in [("fp32", torch.float32), ("fp16", torch.float16), ("bf16", torch.bfloat16)]:
    model = AutoModelForMaskedLM.from_pretrained(MD, torch_dtype=dtype).to(device).eval()
    print(f"\n=== {dtype_name} ===")
    for B in (1, 4, 8, 16, 32, 48):
        x = torch.tensor([ids1] * B, device=device)
        try:
            torch.cuda.reset_peak_memory_stats()
            # warmup
            with torch.no_grad():
                model(input_ids=x[:1])
            torch.cuda.synchronize()
            t0 = time.time()
            with torch.no_grad():
                y = model(input_ids=x).logits
            torch.cuda.synchronize()
            dt = time.time() - t0
            peak = torch.cuda.max_memory_allocated() / 2**30
            del y
            torch.cuda.empty_cache()
            print(f"  B={B:3d}  {dt:7.3f}s  {dt/B:6.3f}s/seq  peak={peak:5.2f} GB")
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            print(f"  B={B:3d}  OOM")
            break
    del model
    torch.cuda.empty_cache()
