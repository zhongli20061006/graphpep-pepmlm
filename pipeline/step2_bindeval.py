#!/usr/bin/env python3
"""Step 2 - score every candidate peptide with BindEvaluator (moPPIt).

BindEvaluator is a *sequence-only* model: it takes the target sequence and a binder
sequence and predicts, for every residue of the TARGET, a sigmoid probability of being
part of the binding site. It does not emit a single peptide-level "binding probability",
so `bind_eval_score` is defined here as:

    binding_sites = { i : p_i >= 0.5 }        (0.5 is the threshold the authors'
                                               own main() hard-codes)
    bind_eval_score = mean(p_i for i in binding_sites),  or 0.0 if that set is empty

That mirrors the authors' own motif-score aggregation (`sum(prediction[motifs])/
len(motifs)`) over the model's own predicted site set. Supporting columns
(`be_n_sites`, `be_max_prob`, `be_mean_prob`) are written to a diagnostics file so the
choice can be audited or re-derived without re-running the model.

Deviations from the reference script (documented):
  * ESM-2 is loaded from the locally converted checkpoint (huggingface.co is unreachable)
  * binders are evaluated in batches; the reference scores one binder per process
  * results are cached by peptide sequence
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import pandas as pd
import torch
from transformers import AutoTokenizer

ROOT = "/home/g/pepmlm_pipeline"
MOPPIT = os.path.join(ROOT, "models", "moPPIt")
ESM2_LOCAL = os.path.join(ROOT, "models", "esm2_t33_650M_UR50D")
CKPT = os.path.join(MOPPIT, "classifier_ckpt", "finetuned_BindEvaluator.ckpt")
CACHE_FILE = os.path.join(ROOT, "cache", "bindeval_cache.json")
ERR_LOG = os.path.join(ROOT, "output", "error.log")

SITE_THRESHOLD = 0.5   # hard-coded in the authors' bindevaluator.py main()


def log_err(msg: str) -> None:
    os.makedirs(os.path.dirname(ERR_LOG), exist_ok=True)
    with open(ERR_LOG, "a") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
    print(f"ERROR: {msg}", file=sys.stderr)


def read_fasta(path: str) -> str:
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line)
    return "".join(parts).upper()


def build_model(device):
    """Rebuild PeptideModel with the authors' published hyper-parameters and load weights."""
    sys.path.insert(0, MOPPIT)
    import torch.nn as nn
    from transformers import EsmModel
    from modules.bindevaluator_modules.models import RepeatedModule3
    from modules.bindevaluator_modules.modules import MultiHeadAttentionSequence, FFN

    class PeptideModel(nn.Module):
        """Exact copy of PeptideModel in bindevaluator.py, with the ESM-2 load made local."""

        def __init__(self, n_layers, d_model, d_hidden, n_head,
                     d_k, d_v, d_inner, dropout=0.2):
            super().__init__()
            self.esm_model = EsmModel.from_pretrained(ESM2_LOCAL)
            for p in self.esm_model.parameters():
                p.requires_grad = False
            self.repeated_module = RepeatedModule3(n_layers, d_model, d_hidden,
                                                   n_head, d_k, d_v, d_inner, dropout=dropout)
            self.final_attention_layer = MultiHeadAttentionSequence(n_head, d_model, d_k, d_v, dropout=dropout)
            self.final_ffn = FFN(d_model, d_inner, dropout=dropout)
            self.output_projection_prot = nn.Linear(d_model, 1)

        def forward(self, binder_tokens, target_tokens):
            peptide_sequence = self.esm_model(**binder_tokens).last_hidden_state
            protein_sequence = self.esm_model(**target_tokens).last_hidden_state
            (prot_enc, sequence_enc, *_rest) = self.repeated_module(peptide_sequence, protein_sequence)
            prot_enc, _ = self.final_attention_layer(prot_enc, sequence_enc, sequence_enc)
            prot_enc = self.final_ffn(prot_enc)
            return self.output_projection_prot(prot_enc)

    # Published hyper-parameters from scripts/bindevaluator.sh
    kwargs = dict(n_layers=8, d_model=128, d_hidden=128, n_head=8,
                  d_k=64, d_v=128, d_inner=64)

    print(f"  torch.load({CKPT}) ...", flush=True)
    ck = torch.load(CKPT, map_location="cpu", weights_only=False)
    print(f"  checkpoint type: {type(ck).__name__}")
    if isinstance(ck, dict) and "state_dict" in ck:
        sd = {k: v for k, v in ck["state_dict"].items()}
        hp = ck.get("hyper_parameters", {}) or {}
        print(f"  checkpoint hparams: {hp}")
    else:
        sd = ck

    model = PeptideModel(**kwargs)
    missing, unexpected = model.load_state_dict(sd, strict=False)
    print(f"  load_state_dict: missing={len(missing)} unexpected={len(unexpected)}")
    if missing:
        print(f"    first missing: {missing[:5]}")
    if unexpected:
        print(f"    first unexpected: {unexpected[:5]}")
    if len(missing) > 0:
        raise RuntimeError(f"checkpoint does not match architecture ({len(missing)} missing)")
    return model.to(device).eval()


def score_batch(model, tok, target_seq, binders, device):
    """Return per-binder (n_sites, mean_site_prob, max_prob, mean_prob, sites, probs)."""
    tgt = tok(target_seq, return_tensors="pt", padding=True, truncation=True, max_length=40000)
    tgt["attention_mask"][0][0] = 0
    tgt["attention_mask"][0][-1] = 0
    tgt_tokens = {k: v.to(device) for k, v in tgt.items()}

    bnd = tok(binders, return_tensors="pt", padding=True, truncation=True, max_length=40000)
    bnd["attention_mask"][:, 0] = 0
    bnd["attention_mask"][:, -1] = 0
    bnd_tokens = {k: v.to(device) for k, v in bnd.items()}

    with torch.no_grad():
        pred = model(bnd_tokens, tgt_tokens).squeeze(-1)   # (B, L_target + 2)
        pred = torch.sigmoid(pred)
    # strip the special-token columns, exactly as the authors do: [0][1:-1]
    probs = pred[:, 1:-1]                                   # (B, L_target)

    out = []
    T = probs.shape[1]
    for row in probs.cpu():
        sites = (row >= SITE_THRESHOLD).nonzero(as_tuple=True)[0]
        n_sites = int(sites.numel())
        score = float(row[sites].mean()) if n_sites > 0 else 0.0
        out.append({
            "n_sites": n_sites,
            "score": score,
            "max_prob": float(row.max()),
            "mean_prob": float(row.mean()),
            "site_fraction": n_sites / T if T else 0.0,
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-fasta", default=os.path.join(ROOT, "data/targets/NCAM1.fasta"))
    ap.add_argument("--peptides", default=os.path.join(ROOT, "output/generated_peptides.csv"))
    ap.add_argument("--out", default=os.path.join(ROOT, "output/generated_peptides.csv"))
    ap.add_argument("--top-ppl", default=os.path.join(ROOT, "output/top20_ppl.txt"))
    ap.add_argument("--top-bindeval", default=os.path.join(ROOT, "output/top20_bindeval.txt"))
    ap.add_argument("--batch-size", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="debug: only score first N")
    args = ap.parse_args()

    t0 = time.time()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    target = read_fasta(args.target_fasta)
    print(f"target length={len(target)}  device={device}  threshold={SITE_THRESHOLD}")

    df = pd.read_csv(args.peptides)
    if args.limit:
        df = df.head(args.limit).copy()
    print(f"candidates={len(df)}")

    print("building BindEvaluator ...", flush=True)
    model = build_model(device)
    tok = AutoTokenizer.from_pretrained(ESM2_LOCAL)

    cache = {}
    if os.path.isfile(CACHE_FILE):
        with open(CACHE_FILE) as fh:
            cache = json.load(fh)
        print(f"  cache entries: {len(cache)}")

    todo = [s for s in dict.fromkeys(df["peptide_sequence"]) if s not in cache]
    print(f"[1/2] scoring {len(todo)} uncached of {df['peptide_sequence'].nunique()} unique", flush=True)

    ts = time.time()
    for i in range(0, len(todo), args.batch_size):
        chunk = todo[i:i + args.batch_size]
        try:
            res = score_batch(model, tok, target, chunk, device)
        except torch.cuda.OutOfMemoryError:
            torch.cuda.empty_cache()
            log_err(f"OOM at batch {i}, retrying singleton for {len(chunk)} peptides")
            res = []
            for s in chunk:
                res.append(score_batch(model, tok, target, [s], device)[0])
                torch.cuda.empty_cache()
        for s, r in zip(chunk, res):
            cache[s] = r
        done = min(i + args.batch_size, len(todo))
        if (i // args.batch_size) % 10 == 0 or done == len(todo):
            print(f"      {done}/{len(todo)}  {time.time()-ts:.1f}s", flush=True)
        os.makedirs(os.path.dirname(CACHE_FILE), exist_ok=True)
        with open(CACHE_FILE, "w") as fh:
            json.dump(cache, fh)
    print(f"      scoring done in {time.time()-ts:.1f}s")

    df["bind_eval_score"] = [cache[s]["score"] for s in df["peptide_sequence"]]

    diag = df[["peptide_id", "peptide_sequence"]].copy()
    for col, key in [("be_n_sites", "n_sites"), ("be_max_prob", "max_prob"),
                     ("be_mean_prob", "mean_prob"), ("be_site_fraction", "site_fraction")]:
        diag[col] = [cache[s][key] for s in df["peptide_sequence"]]
    diag_path = os.path.join(ROOT, "output", "bindeval_diagnostics.csv")
    diag.to_csv(diag_path, index=False)

    # ---- choose the ranking column, then persist the scored table ----
    # The task specifies bind_eval_score. If it is constant across the pool (which happens
    # when no residue ever crosses the site threshold, or when every peptide hits the same
    # site set) it cannot rank anything, so we fall back to the max site probability and
    # record that choice in the table and the report.
    be_unique = df["bind_eval_score"].nunique(dropna=True)
    rank_col = "bind_eval_score"
    if be_unique <= 1:
        rank_col = "be_max_prob"
        msg = (f"bind_eval_score has only {be_unique} distinct value(s) across {len(df)} "
               f"candidates -> cannot discriminate; ranking the BindEvaluator arm by "
               f"be_max_prob instead (documented in the report).")
        log_err(msg)
        print(f"  WARNING: {msg}")
        rank_series = diag["be_max_prob"].to_numpy()
    else:
        rank_series = df["bind_eval_score"].to_numpy()

    print(f"  ranking column for the BindEvaluator arm: {rank_col}")
    df = df.assign(bind_eval_rank_col=rank_col)
    df.to_csv(args.out, index=False)

    # ---- Top-20 lists ----
    rank_df = df.assign(_rank=rank_series)
    top_ppl = df.nsmallest(20, "ppl_score")
    top_be = rank_df.nlargest(20, "_rank").drop(columns=["_rank"])

    def write_list(sub, path):
        with open(path, "w") as fh:
            fh.write("peptide_id\tpeptide_sequence\tppl_score\tbind_eval_score\n")
            for _, r in sub.iterrows():
                fh.write(f"{r['peptide_id']}\t{r['peptide_sequence']}\t"
                         f"{r['ppl_score']:.6f}\t{r['bind_eval_score']:.6f}\n")

    write_list(top_ppl, args.top_ppl)
    write_list(top_be, args.top_bindeval)

    # combined, de-duplicated list for the structure stage (one row per unique peptide)
    comb = pd.concat([
        top_ppl.assign(source="ppl"),
        top_be.assign(source="bindeval"),
    ], ignore_index=True)
    comb = comb.drop_duplicates(subset=["peptide_sequence"], keep="first")
    comb_out = os.path.join(ROOT, "output", "eval_lists.csv")
    comb[["peptide_id", "peptide_sequence", "ppl_score", "bind_eval_score", "source"]] \
        .to_csv(comb_out, index=False)
    print(f"  {comb_out}  ({len(comb)} unique peptides for the structure stage)")

    dt = time.time() - t0
    print(f"\nSTEP2 DONE  time={dt:.1f}s")
    print(f"  {args.out}   (added bind_eval_score)")
    print(f"  {diag_path}   (diagnostics)")
    print(f"  {args.top_ppl}   top20 by ppl ascending")
    print(f"  {args.top_bindeval}   top20 by bind_eval_score descending")
    print(f"  bind_eval_score: min={df['bind_eval_score'].min():.4f} "
          f"median={df['bind_eval_score'].median():.4f} max={df['bind_eval_score'].max():.4f}")
    n_zero = int((df["bind_eval_score"] == 0).sum())
    print(f"  peptides with an EMPTY predicted binding site (score=0): {n_zero}/{len(df)}")
    print(f"\ntop5 by bind_eval_score:\n{top_be.head(5).to_string(index=False)}")
    print(f"\ntop5 by ppl:\n{top_ppl.head(5).to_string(index=False)}")
    overlap = set(top_ppl["peptide_sequence"]) & set(top_be["peptide_sequence"])
    print(f"\noverlap between the two Top-20 lists: {len(overlap)} sequences")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
