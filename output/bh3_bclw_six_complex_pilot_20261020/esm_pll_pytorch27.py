"""PyTorch 2.6+ safe-loader compatibility entry point for the frozen ESM PLL."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import esm_pll


def load_local_checkpoint(checkpoint: Path):
    """Keep weights_only safe loading; allow the one global present in this checkpoint."""
    import argparse as argparse_module
    import esm
    import torch

    with torch.serialization.safe_globals([argparse_module.Namespace]):
        return esm.pretrained.load_model_and_alphabet_local(str(checkpoint))


def run(input_fasta: Path, output_csv: Path, checkpoint: Path, device_name: str) -> None:
    import torch

    device = "cuda" if device_name == "auto" and torch.cuda.is_available() else device_name
    if device == "auto":
        device = "cpu"
    model, alphabet = load_local_checkpoint(checkpoint)
    model.eval().to(device)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sequence_id", "sequence", "sum_masked_log_likelihood", "mean_masked_log_probability", "position_log_probabilities"),
        )
        writer.writeheader()
        for identifier, sequence in esm_pll.read_fasta(input_fasta):
            position_scores = esm_pll.score_sequence(model, alphabet, sequence, device)
            writer.writerow({
                "sequence_id": identifier,
                "sequence": sequence,
                "sum_masked_log_likelihood": f"{sum(position_scores):.12g}",
                "mean_masked_log_probability": f"{sum(position_scores) / len(position_scores):.12g}",
                "position_log_probabilities": ";".join(f"{value:.12g}" for value in position_scores),
            })


def main() -> None:
    parser = argparse.ArgumentParser(description="ESM-2 masked pseudo-log-likelihood with PyTorch 2.6+ safe-loader compatibility")
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    args = parser.parse_args()
    run(args.input_fasta, args.output_csv, args.checkpoint, args.device)


if __name__ == "__main__":
    main()
