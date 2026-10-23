#!/usr/bin/env python3
"""Step 1 - generate candidate peptides with PepMLM-650M.

Faithful to the authors' official implementation (HuggingFace Space
`TianlaiChen/PepMLM` app.py):

  * input  = target_sequence + "<mask>" * peptide_length   (single forward pass)
  * sample = top-k over the logits at the mask positions, then one categorical draw
  * PPL    = "pseudo perplexity": for every binder position, re-mask that position
             alone, take the cross-entropy of the TRUE token, average, then exp()

Deviations (documented, deliberate):
  * batched sampling + chunked pseudo-PPL scoring (the reference loops one peptide
    at a time, which is far too slow for 1000 candidates on this GPU)
  * an on-disk cache keyed by peptide sequence so re-runs never recompute PPL

Output: output/generated_peptides.csv  (peptide_id, peptide_sequence, ppl_score)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.distributions.categorical import Categorical
from transformers import AutoModelForMaskedLM, AutoTokenizer

ROOT = "/home/g/pepmlm_pipeline"
MODEL_DIR = os.path.join(ROOT, "models", "PepMLM-650M")
CACHE_DIR = os.path.join(ROOT, "cache")
CACHE_FILE = os.path.join(CACHE_DIR, "ppl_cache.json")
ERR_LOG = os.path.join(ROOT, "output", "error.log")


def log_err(msg: str) -> None:
    os.makedirs(os.path.dirname(ERR_LOG), exist_ok=True)
    with open(ERR_LOG, "a") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
    print(f"ERROR: {msg}", file=sys.stderr)


def read_fasta(path: str) -> tuple[str, str]:
    """Return (header, sequence) from a single-record FASTA."""
    header, parts = "", []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    break
                header = line[1:]
            else:
                parts.append(line)
    seq = "".join(parts).upper()
    if not seq:
        raise ValueError(f"no sequence found in {path}")
    return header, seq


def load_cache() -> dict:
    if os.path.isfile(CACHE_FILE):
        try:
            with open(CACHE_FILE) as fh:
                return json.load(fh)
        except Exception as exc:  # noqa: BLE001
            log_err(f"ppl cache unreadable, starting fresh: {exc}")
    return {}


def save_cache(cache: dict) -> None:
    os.makedirs(CACHE_DIR, exist_ok=True)
    tmp = CACHE_FILE + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(cache, fh)
    os.replace(tmp, CACHE_FILE)


def sample_peptides(model, tok, target_seq, length, top_k, n, batch_size, device, seed):
    """One forward pass per batch; top-k categorical sampling at mask positions."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    seqs: list[str] = []
    input_seq = target_seq + "<mask>" * length
    enc = tok(input_seq, return_tensors="pt")
    n_tok = int(enc["input_ids"].shape[1])

    made = 0
    while made < n:
        b = min(batch_size, n - made)
        ids = enc["input_ids"].repeat(b, 1).to(device)
        attn = enc["attention_mask"].repeat(b, 1).to(device)
        with torch.no_grad():
            logits = model(input_ids=ids, attention_mask=attn).logits
        mask_pos = (ids[0] == tok.mask_token_id).nonzero(as_tuple=True)[0]
        if mask_pos.numel() != length:
            raise RuntimeError(f"expected {length} mask positions, found {mask_pos.numel()}")
        # (b, length, vocab)
        at_masks = logits[:, mask_pos, :]
        topk_logits, topk_idx = at_masks.topk(top_k, dim=-1)
        probs = F.softmax(topk_logits, dim=-1)
        draws = Categorical(probs.reshape(-1, top_k)).sample().reshape(b, length, 1)
        chosen = topk_idx.gather(-1, draws).squeeze(-1)
        for row in chosen.cpu().tolist():
            s = "".join(tok.convert_ids_to_tokens(i) for i in row)
            seqs.append(s)
        made += b
        del logits, at_masks, topk_logits, topk_idx, probs, draws, chosen
        if device == "cuda":
            torch.cuda.empty_cache()
    return seqs, n_tok


def pseudo_ppl_chunk(model, tok, target_ids, peptides, device, n_pos, row_chunk):
    """Authors' pseudo-perplexity, computed with one forward pass per mask position.

    For a chunk of B peptides we build B masked variants per binder position j: the full
    sequence is reconstructed as target_ids + peptide_residues + [<eos>] (matching
    `tokenize(target + peptide)` exactly), with binder position j replaced by <mask>.

    Reference behaviour reproduced exactly: mask ONE binder position, take the
    cross-entropy of the TRUE token at that position, average over the n_pos positions,
    then exp().
    """
    B = len(peptides)
    pep_ids = [tok(p, add_special_tokens=False)["input_ids"] for p in peptides]
    L = target_ids.shape[0]                    # <cls> + target residues
    eos_id = tok.eos_token_id
    if any(len(p) != n_pos for p in pep_ids):
        raise ValueError(f"all peptides must tokenise to {n_pos} residues: "
                         f"{[len(p) for p in pep_ids]}")

    losses = torch.empty(B, n_pos)
    for j in range(n_pos):          # outer loop = one masked column per pass
        rows = []
        for b in range(B):
            row = list(pep_ids[b])
            row[j] = tok.mask_token_id
            rows.append(target_ids.tolist() + row + [eos_id])
        true_j = torch.tensor([pep_ids[b][j] for b in range(B)], device=device)
        row_losses = torch.empty(B, device=device)
        for s in range(0, B, row_chunk):
            sub = torch.tensor(rows[s:s + row_chunk], device=device)
            with torch.no_grad():
                logits = model(input_ids=sub).logits
            # masked column j of the binder sits at index L + j
            logp = F.log_softmax(logits[:, L + j, :].float(), dim=-1)
            row_losses[s:s + sub.shape[0]] = -logp.gather(
                1, true_j[s:s + sub.shape[0]].unsqueeze(1)).squeeze(1)
            del logits, logp
            if device == "cuda":
                torch.cuda.empty_cache()
        losses[:, j] = row_losses.cpu()

    return torch.exp(losses.mean(dim=1)).tolist()


def pseudo_ppl_reference(model, tok, target_seq, peptide, device):
    """Verbatim port of the authors' per-peptide loop, used to validate the fast path."""
    n_pos = len(peptide)
    tensor_input = tok.encode(target_seq + peptide, return_tensors="pt").to(device)
    total = 0.0
    for i in range(-n_pos - 1, -1):
        masked_input = tensor_input.clone()
        masked_input[0, i] = tok.mask_token_id
        labels = torch.full(tensor_input.shape, -100).to(device)
        labels[0, i] = tensor_input[0, i]
        with torch.no_grad():
            total += model(masked_input, labels=labels).loss.item()
    return float(np.exp(total / n_pos))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-fasta", default=os.path.join(ROOT, "data/targets/NCAM1.fasta"))
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--length", type=int, default=10)
    ap.add_argument("--top-k", type=int, default=3)
    ap.add_argument("--sample-batch", type=int, default=50)
    ap.add_argument("--ppl-chunk", type=int, default=32, help="peptides per pseudo-PPL chunk")
    ap.add_argument("--ppl-row-chunk", type=int, default=8, help="rows per forward pass")
    ap.add_argument("--verify-ppl", type=int, default=3,
                    help="compare N peptides against the authors' reference loop")
    ap.add_argument("--dtype", choices=["fp32", "fp16", "bf16"], default="fp32",
                    help="inference precision; fp16/bf16 are faster but change PPL slightly")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default=os.path.join(ROOT, "output/generated_peptides.csv"))
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    t0 = time.time()
    device = "cuda" if torch.cuda.is_available() else "cpu"

    header, target = read_fasta(args.target_fasta)
    print(f"target: {header[:60]}...  length={len(target)}")
    print(f"device={device}  n={args.n}  peptide_length={args.length}  top_k={args.top_k}")

    print("loading PepMLM-650M ...", flush=True)
    tok = AutoTokenizer.from_pretrained(MODEL_DIR)
    dtype = {"fp32": torch.float32, "fp16": torch.float16, "bf16": torch.bfloat16}[args.dtype]
    model = AutoModelForMaskedLM.from_pretrained(MODEL_DIR, torch_dtype=dtype).to(device).eval()
    n_param = sum(p.numel() for p in model.parameters())
    print(f"  params={n_param/1e6:.1f}M  vocab={tok.vocab_size}  dtype={args.dtype}")
    if args.dtype != "fp32":
        print("  NOTE: reduced precision -> PPL is computed in fp32 log-softmax but activations "
              "are lower precision; the verification below quantifies the deviation.")

    seq_tok_len = len(tok(target + "<mask>" * args.length)["input_ids"])
    print(f"  tokenised input length = {seq_tok_len} (ESM-2 max_position_embeddings="
          f"{model.config.max_position_embeddings})")
    if seq_tok_len > model.config.max_position_embeddings:
        log_err(f"input length {seq_tok_len} exceeds ESM-2 position limit; will be truncated")

    # target tokenised once; peptide tokens are appended directly. We must replicate
    # `tokenize(target + peptide)` exactly:
    #   [<cls>] + target_residues + peptide_residues + [<eos>]
    # NOTE: tok(peptide) would add its own <cls>/<eos>, so special tokens are suppressed
    # and the single trailing <eos> is appended by hand. (Verified against
    # tok(target+peptide) token-for-token in scripts/check_tok_align.py.)
    target_full = tok(target)["input_ids"]
    assert target_full[0] == tok.cls_token_id and target_full[-1] == tok.eos_token_id
    target_ids = torch.tensor(target_full[:-1] )          # <cls> + target residues (no <eos>)
    eos_id = tok.eos_token_id

    ref_ids = tok(target + "A" * args.length)["input_ids"]
    probe_ids = target_ids.tolist() + [tok.convert_tokens_to_ids("A")] * args.length + [eos_id]
    assert probe_ids == ref_ids, "token assembly does not match tok(target+peptide)"
    print(f"  target tokens={len(target_ids)} (+1 eos)  binder tokens={args.length}  "
          f"total={len(ref_ids)}  [token assembly verified]")

    # ---- sample ----
    print(f"\n[1/2] sampling {args.n} peptides ...", flush=True)
    ts = time.time()
    peptides, _ = sample_peptides(model, tok, target, args.length, args.top_k,
                                  args.n, args.sample_batch, device, args.seed)
    print(f"      sampled {len(peptides)} in {time.time()-ts:.1f}s")

    # ---- pseudo-perplexity with cache ----
    cache = load_cache()
    uniq_missing = sorted({p for p in peptides if p not in cache})
    print(f"[2/2] pseudo-PPL: {len(set(peptides))} unique, {len(uniq_missing)} uncached", flush=True)

    # validate the fast path against the authors' literal loop on a few peptides
    if args.verify_ppl > 0 and uniq_missing:
        check = uniq_missing[:args.verify_ppl]
        fast = pseudo_ppl_chunk(model, tok, target_ids, check, device, args.length, args.ppl_row_chunk)
        ref = [pseudo_ppl_reference(model, tok, target, p, device) for p in check]
        print("      fast-path vs authors' reference loop:")
        worst = 0.0
        for s, a, b in zip(check, fast, ref):
            rel = abs(a - b) / max(abs(b), 1e-9)
            worst = max(worst, rel)
            print(f"        {s}  fast={a:.6f}  ref={b:.6f}  rel={rel:.2e}")
        # fp32 must match the authors' loop to numerical precision; reduced precision
        # carries an inherent ~1e-4..1e-3 deviation, which we bound and record.
        tol = 1e-4 if args.dtype == "fp32" else 5e-3
        if worst > tol:
            log_err(f"pseudo-PPL fast path disagrees with reference "
                    f"(dtype={args.dtype}, worst rel={worst:.3e} > tol={tol:.0e})")
            return 3
        print(f"      PASS (worst relative difference {worst:.2e}, tolerance {tol:.0e}, "
              f"dtype={args.dtype})")
    ts = time.time()
    for i in range(0, len(uniq_missing), args.ppl_chunk):
        chunk = uniq_missing[i:i + args.ppl_chunk]
        vals = pseudo_ppl_chunk(model, tok, target_ids, chunk, device,
                                args.length, args.ppl_row_chunk)
        for s, v in zip(chunk, vals):
            cache[s] = float(v)
        done = min(i + args.ppl_chunk, len(uniq_missing))
        print(f"      PPL {done}/{len(uniq_missing)}  ({time.time()-ts:.1f}s)", flush=True)
        save_cache(cache)
    print(f"      PPL stage done in {time.time()-ts:.1f}s")

    df = pd.DataFrame({
        "peptide_id": [f"PEP{i:05d}" for i in range(len(peptides))],
        "peptide_sequence": peptides,
        "ppl_score": [cache[p] for p in peptides],
    })
    df.to_csv(args.out, index=False)

    dt = time.time() - t0
    print(f"\nSTEP1 DONE  time={dt:.1f}s  file={args.out}")
    print(f"  rows={len(df)}  unique_sequences={df['peptide_sequence'].nunique()}")
    print(f"  ppl: min={df['ppl_score'].min():.3f} median={df['ppl_score'].median():.3f} "
          f"max={df['ppl_score'].max():.3f}")
    print(f"  top5 by PPL (lowest):\n{df.nsmallest(5, 'ppl_score').to_string(index=False)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
