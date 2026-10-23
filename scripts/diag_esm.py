#!/usr/bin/env python3
"""Diagnose why the converted model predicts the wrong masked token."""
import torch
from transformers import AutoTokenizer, EsmForMaskedLM, EsmModel

LOCAL = "/home/g/pepmlm_pipeline/models/esm2_t33_650M_UR50D"
FAIR = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"

tok = AutoTokenizer.from_pretrained(LOCAL)
print("=== tokenizer ===")
print("vocab size:", tok.vocab_size)
print("mask id:", tok.mask_token_id, "pad:", tok.pad_token_id,
      "cls:", tok.cls_token_id, "eos:", tok.eos_token_id)
print("first 10 tokens:", tok.convert_ids_to_tokens(list(range(10))))
seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYG"
ids = tok(seq)["input_ids"]
print("encoded len:", len(ids), "->", tok.convert_ids_to_tokens(ids[:12]))

print("\n=== local model ===")
mdl = EsmForMaskedLM.from_pretrained(LOCAL).eval()
print("lm_head.weight sum:", float(torch.load(LOCAL + "/pytorch_model.bin",
      map_location="cpu", weights_only=True, mmap=True)["lm_head.decoder.weight"].float().sum()))
print("word_emb sum:", float(mdl.esm.embeddings.word_embeddings.weight.float().sum()))

masked = seq[:10] + "<mask>" + seq[11:]
m = tok(masked, return_tensors="pt")
print("masked ids:", m["input_ids"].tolist())
pos = (m["input_ids"] == tok.mask_token_id).nonzero()[0, 1].item()
print("mask pos:", pos)
with torch.no_grad():
    lo = mdl(**m).logits[0, pos]
top = torch.topk(lo, 8)
print("top8:", [(tok.convert_ids_to_tokens([i.item()])[0], round(v.item(), 2))
                for v, i in zip(top.values, top.indices)])
print("logit of true token:", round(lo[ids[pos]].item(), 3),
      "rank:", int((lo > lo[ids[pos]]).sum().item()))

print("\n=== compare against a genuine transformers ESM-2 (random init would be ~uniform) ===")
# Sanity: what does the *same* architecture with random weights predict? (baseline for comparison)
torch.manual_seed(0)
import json
from transformers import EsmConfig
with open(LOCAL + "/config.json") as fh:
    cfg = EsmConfig(**json.load(fh))
rnd = EsmForMaskedLM(cfg).eval()
with torch.no_grad():
    lo2 = rnd(**m).logits[0, pos]
t2 = torch.topk(lo2, 5)
print("random-init top5:", [(tok.convert_ids_to_tokens([i.item()])[0], round(v.item(), 2))
                            for v, i in zip(t2.values, t2.indices)])

print("\n=== per-layer energy check (are weights actually loaded?) ===")
sd = torch.load(LOCAL + "/pytorch_model.bin", map_location="cpu", weights_only=True, mmap=True)
for k in ["esm.embeddings.word_embeddings.weight",
          "esm.encoder.layer.0.attention.self.query.weight",
          "esm.encoder.emb_layer_norm_after.weight",
          "lm_head.decoder.weight"]:
    t = sd[k].float()
    print(f"  {k:55s} mean={t.mean():+.4f} std={t.std():.4f}")

# Compare to the source checkpoint directly
print("\n=== source fair-esm values for the same tensors ===")
obj = torch.load(FAIR, map_location="cpu", weights_only=False, mmap=True)
s = obj["model"]
for k in ["encoder.sentence_encoder.embed_tokens.weight",
          "encoder.sentence_encoder.layers.0.self_attn.q_proj.weight",
          "encoder.sentence_encoder.emb_layer_norm_after.weight",
          "encoder.lm_head.weight"]:
    t = s[k].float()
    print(f"  {k:60s} mean={t.mean():+.4f} std={t.std():.4f}")
