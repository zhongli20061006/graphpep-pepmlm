#!/usr/bin/env python3
"""Verify token alignment between tokenize(target+pep) and tokenize(target)+tokenize(pep),
and check that the binder occupies exactly the last L tokens with no boundary merging.
"""
from transformers import AutoTokenizer

MD = "/home/g/pepmlm_pipeline/models/PepMLM-650M"
tok = AutoTokenizer.from_pretrained(MD)

tgt = ("MLQTKDLIWTLFFLGTAVSLQVDIVPSQGEISVGESKFFLCQVAGDAKDKDISWFSPNGEKLTPNQQRISV"
       "VWNDDSSTSKLQISNVTDDDGGIYKCIVTGEDKDQQQ")  # arbitrary short fragment for clarity
pep = "KIKKLLARIT"
L = len(pep)

full = tok(tgt + pep)["input_ids"]
t_only = tok(tgt)["input_ids"]
p_only = tok(pep)["input_ids"]
print("len(full)      =", len(full))
print("len(target)    =", len(t_only))
print("len(pep)       =", len(p_only))
print("target tail    =", t_only[-4:], tok.convert_ids_to_tokens(t_only[-4:]))
print("pep tokens     =", p_only, tok.convert_ids_to_tokens(p_only))
print("full tail(L+3) =", full[-(L + 3):], tok.convert_ids_to_tokens(full[-(L + 3):]))
print()
print("full[-(L+1):-1] == pep tokens ?", full[-(L + 1):-1] == p_only)

# the real case uses the full 858-residue NCAM1
import os
fasta = "/home/g/pepmlm_pipeline/data/targets/NCAM1.fasta"
seq = "".join(l.strip() for l in open(fasta) if not l.startswith(">")).upper()
full2 = tok(seq + pep)["input_ids"]
t2 = tok(seq)["input_ids"]
print("\n=== real NCAM1 (858 aa) ===")
print("len(target tokens) =", len(t2))
print("len(full tokens)   =", len(full2))
print("expected len       =", len(t2) + L)
print("binder region      =", full2[-(L + 1):-1], "==", p_only, "?", full2[-(L + 1):-1] == p_only)

# masked form
mask_full = tok(seq + "<mask>" * L)["input_ids"]
print("\n=== masked form ===")
print("len(masked)        =", len(mask_full))
print("mask ids           =", mask_full[-(L + 1):-1])
print("mask_token_id      =", tok.mask_token_id)
print("all mask?          =", all(i == tok.mask_token_id for i in mask_full[-(L + 1):-1]))
print("last token (eos?)  =", mask_full[-1], tok.convert_ids_to_tokens([mask_full[-1]]))
