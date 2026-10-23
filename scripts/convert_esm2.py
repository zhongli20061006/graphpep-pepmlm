#!/usr/bin/env python3
"""Convert the local fair-esm ESM-2 650M checkpoint into `transformers` format.

Why: BindEvaluator calls `EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D")`,
which needs huggingface.co — unreachable on this machine. The identical ESM-2 650M weights
already exist locally in fair-esm format, so we re-serialise them instead of downloading
2.6 GB again.

Mapping was derived empirically: instantiate `EsmForMaskedLM` from the bundled config, dump
its `state_dict()` keys, and diff against the fair-esm checkpoint keys. Rules:
  encoder.sentence_encoder.layers.N.<suffix>  -> esm.encoder.layer.N.<mapped suffix>
  encoder.sentence_encoder.embed_tokens.*     -> esm.embeddings.word_embeddings.*
  encoder.sentence_encoder.<rest>             -> esm.<rest>          (drop sentence_encoder)
  encoder.lm_head.<rest>                      -> lm_head.<rest>
  encoder.emb_layer_norm_after.<rest>         -> esm.encoder.emb_layer_norm_after.<rest>

Output: /home/g/pepmlm_pipeline/models/esm2_t33_650M_UR50D/
"""
from __future__ import annotations

import gc
import json
import os
import shutil
import sys

import torch

FAIR_ESM_CKPT = "/home/g/GraphPep_v1.3/models/esm2_t33_650M_UR50D.pt"
OUT_DIR = "/home/g/pepmlm_pipeline/models/esm2_t33_650M_UR50D"
PEPMLM_DIR = "/home/g/pepmlm_pipeline/models/PepMLM-650M"

# Per-layer suffix mapping: fair-esm -> transformers
LAYER_SUFFIX = {
    "self_attn.k_proj": "attention.self.key",
    "self_attn.v_proj": "attention.self.value",
    "self_attn.q_proj": "attention.self.query",
    "self_attn.out_proj": "attention.output.dense",
    "self_attn.rot_emb": "attention.self.rotary_embeddings",
    "self_attn_layer_norm": "attention.LayerNorm",
    "fc1": "intermediate.dense",
    "fc2": "output.dense",
    "final_layer_norm": "LayerNorm",
}

SE_SENT = "encoder.sentence_encoder."


def map_key(k: str) -> str | None:
    """Return the transformers key for a fair-esm key, or None if it must be dropped."""
    if k.startswith(SE_SENT + "layers."):
        rest = k[len(SE_SENT + "layers."):]
        idx, _, suffix = rest.partition(".")
        for fsrc, fdst in LAYER_SUFFIX.items():
            if suffix.startswith(fsrc):
                return f"esm.encoder.layer.{idx}.{fdst}{suffix[len(fsrc):]}"
        raise KeyError(f"unmapped layer suffix: {k}")
    if k == SE_SENT + "embed_tokens.weight":
        return "esm.embeddings.word_embeddings.weight"
    if k.startswith(SE_SENT + "emb_layer_norm_after."):
        # transformers keeps this one under esm.encoder (not esm.)
        return "esm.encoder." + k[len(SE_SENT):]
    if k.startswith(SE_SENT):
        return "esm." + k[len(SE_SENT):]
    if k.startswith("encoder.lm_head."):
        rest = k[len("encoder.lm_head."):]
        if rest == "weight":
            # fair-esm ties the MLM decoder to the token table; transformers exposes
            # `lm_head.weight` (ignored at load because lm_head.decoder is the tied one).
            return "lm_head.weight"
        return "lm_head." + rest
    if k.startswith("encoder.emb_layer_norm_after."):
        return "esm.encoder." + k[len("encoder."):]
    if "contact_head" in k:
        return "esm.contact_head." + k.split("contact_head.")[1]
    return None


def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"[1/4] mmap-loading {FAIR_ESM_CKPT}", flush=True)
    obj = torch.load(FAIR_ESM_CKPT, map_location="cpu", weights_only=False, mmap=True)
    src = obj["model"]
    print(f"      {len(src)} source keys", flush=True)

    new_sd: dict[str, torch.Tensor] = {}
    dropped: list[str] = []
    for k, v in src.items():
        nk = map_key(k)
        if nk is None:
            dropped.append(k)
            continue
        new_sd[nk] = v
    if dropped:
        print(f"      dropped {len(dropped)} fair-esm-only keys: {dropped[:6]}")

    # `lm_head.decoder` is the tied MLM decoder in transformers; expose it explicitly so the
    # mapping covers every key of EsmForMaskedLM.state_dict().
    if "esm.embeddings.word_embeddings.weight" in new_sd:
        new_sd["lm_head.decoder.weight"] = new_sd["esm.embeddings.word_embeddings.weight"]

    print(f"      mapped -> {len(new_sd)} keys", flush=True)

    # ---- strict verification against the real transformers model ----
    print("[2/4] verifying against transformers state_dict ...", flush=True)
    from transformers import EsmConfig, EsmForMaskedLM

    with open(os.path.join(PEPMLM_DIR, "config.json")) as fh:
        cfg = EsmConfig(**json.load(fh))
    ref = EsmForMaskedLM(cfg).state_dict()

    missing = [k for k in ref if k not in new_sd]
    # `lm_head.weight` is not part of EsmForMaskedLM.state_dict (it exposes lm_head.decoder),
    # so it is expected to be "extra"; everything else extra is suspicious.
    extra = [k for k in new_sd if k not in ref and k != "lm_head.weight"]
    bad = [
        (k, tuple(new_sd[k].shape), tuple(ref[k].shape))
        for k in ref
        if k in new_sd and tuple(new_sd[k].shape) != tuple(ref[k].shape)
    ]
    print(f"      ref={len(ref)} ours={len(new_sd)} missing={len(missing)} extra={len(extra)} shape_bad={len(bad)}")
    for k in missing[:10]:
        print(f"        MISSING {k} {tuple(ref[k].shape)}")
    for k in extra[:10]:
        print(f"        EXTRA   {k} {tuple(new_sd[k].shape)}")
    for k, a, b in bad[:10]:
        print(f"        SHAPE   {k} ours={a} ref={b}")

    # contact_head is optional (GraphPep/BindEvaluator never use it).
    hard_missing = [k for k in missing if "contact_head" not in k]
    if hard_missing or bad:
        raise SystemExit(f"FATAL: {len(hard_missing)} missing, {len(bad)} shape mismatches")

    # `lm_head.decoder.weight` must equal the token embedding table (weight tying).
    if not torch.equal(new_sd["lm_head.decoder.weight"], new_sd["esm.embeddings.word_embeddings.weight"]):
        raise SystemExit("FATAL: lm_head.decoder.weight is not tied to word_embeddings")

    print("[3/4] saving ...", flush=True)
    torch.save(new_sd, os.path.join(OUT_DIR, "pytorch_model.bin"))
    del new_sd, ref, src, obj
    gc.collect()

    shutil.copy2(os.path.join(PEPMLM_DIR, "config.json"), os.path.join(OUT_DIR, "config.json"))
    for fn in ("vocab.txt", "tokenizer_config.json", "special_tokens_map.json"):
        p = os.path.join(PEPMLM_DIR, fn)
        if os.path.isfile(p):
            shutil.copy2(p, os.path.join(OUT_DIR, fn))

    print("[4/4] round-trip load + functional check ...", flush=True)
    from transformers import AutoTokenizer, EsmForMaskedLM, EsmModel

    tok = AutoTokenizer.from_pretrained(OUT_DIR)
    mdl = EsmForMaskedLM.from_pretrained(OUT_DIR)
    mdl.eval()

    # Weight-tie invariant: the tied decoder must equal the token embedding table.
    assert torch.equal(
        mdl.lm_head.decoder.weight, mdl.esm.embeddings.word_embeddings.weight
    ), "lm_head.decoder.weight is not tied to word_embeddings"

    # Functional check against the original fair-esm implementation: the converted HF
    # model must reproduce fair-esm logits on identical token ids. (A single masked-token
    # argmax is NOT a valid check — ESM-2 650M itself can predict a different residue.)
    seq = "MKTAYIAKQRQISFVKSHFSRQDILDLWIYHTQGYFPDWQNYG"
    masked_seq = seq[:10] + "<mask>" + seq[11:]
    ids = tok(seq, return_tensors="pt")
    m = tok(masked_seq, return_tensors="pt")
    with torch.no_grad():
        h_logits = mdl(**m).logits[0]

    import esm as fair_esm_mod

    _orig = fair_esm_mod.pretrained.torch.load
    fair_esm_mod.pretrained.torch.load = lambda *a, **kw: _orig(*a, **{**kw, "weights_only": False})
    fair_model, alphabet = fair_esm_mod.pretrained.load_model_and_alphabet_local(FAIR_ESM_CKPT)
    fair_esm_mod.pretrained.torch.load = _orig
    fair_model.eval()
    f_tokens = torch.tensor([m["input_ids"][0].tolist()])
    with torch.no_grad():
        f_logits = fair_model(f_tokens)["logits"][0]

    max_diff = (h_logits - f_logits).abs().max().item()
    print(f"      HF vs fair-esm max |logit diff| = {max_diff:.3e}")

    em = EsmModel.from_pretrained(OUT_DIR)
    em.eval()
    with torch.no_grad():
        hs = em(**ids)
    print(f"      EsmModel.last_hidden_state={tuple(hs.last_hidden_state.shape)}")

    if max_diff > 1e-3:
        print("FATAL: converted model does not match fair-esm", file=sys.stderr)
        return 2
    print(f"DONE -> {OUT_DIR}  (numerically equivalent to fair-esm)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
