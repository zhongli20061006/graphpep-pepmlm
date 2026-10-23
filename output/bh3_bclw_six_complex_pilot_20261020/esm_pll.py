from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Sequence

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def validate_sequence(sequence: str) -> str:
    normalized = sequence.strip().upper()
    if not normalized:
        raise ValueError("sequence must not be empty")
    invalid = sorted(set(normalized) - STANDARD_AMINO_ACIDS)
    if invalid:
        raise ValueError(f"non-standard residue(s): {''.join(invalid)}")
    return normalized


def masked_sequences(sequence: str) -> tuple[str, ...]:
    normalized = validate_sequence(sequence)
    return tuple(normalized[:i] + "<mask>" + normalized[i + 1 :] for i in range(len(normalized)))


def sum_masked_log_probs(
    logits: Sequence[Sequence[float]], target_token_ids: Sequence[int]
) -> float:
    if len(logits) != len(target_token_ids):
        raise ValueError("one logit row is required for each target token")
    total = 0.0
    for row, target_id in zip(logits, target_token_ids):
        if not row:
            raise ValueError("logit rows must not be empty")
        if target_id < 0 or target_id >= len(row):
            raise ValueError("target token index is outside the logit row")
        values = [float(value) for value in row]
        if not all(math.isfinite(value) for value in values):
            raise ValueError("logits must be finite")
        peak = max(values)
        log_normalizer = peak + math.log(sum(math.exp(value - peak) for value in values))
        total += values[target_id] - log_normalizer
    return total


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_id: str | None = None
    current_sequence: list[str] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id is not None:
                records.append((current_id, validate_sequence("".join(current_sequence))))
            current_id = line[1:].strip().split()[0] if line[1:].strip() else ""
            if not current_id:
                raise ValueError(f"empty FASTA identifier on line {line_number}")
            current_sequence = []
        elif current_id is None:
            raise ValueError(f"sequence data before FASTA header on line {line_number}")
        else:
            current_sequence.append(line)
    if current_id is not None:
        records.append((current_id, validate_sequence("".join(current_sequence))))
    if not records:
        raise ValueError("FASTA has no records")
    identifiers = [identifier for identifier, _ in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("FASTA identifiers must be unique")
    return records


def score_sequence(model: object, alphabet: object, sequence: str, device: str) -> list[float]:
    import torch
    import torch.nn.functional as functional

    sequence = validate_sequence(sequence)
    masked = masked_sequences(sequence)
    batch_converter = alphabet.get_batch_converter()
    scores: list[float] = []
    with torch.no_grad():
        for position, masked_sequence in enumerate(masked):
            _, _, tokens = batch_converter([(str(position), masked_sequence)])
            logits = model(tokens.to(device))["logits"][0, position + 1]
            target_id = alphabet.get_idx(sequence[position])
            scores.append(float(functional.log_softmax(logits, dim=-1)[target_id].item()))
    return scores


def run(input_fasta: Path, output_csv: Path, checkpoint: Path, device_name: str) -> None:
    import esm
    import torch

    device = "cuda" if device_name == "auto" and torch.cuda.is_available() else device_name
    if device == "auto":
        device = "cpu"
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(str(checkpoint))
    model.eval().to(device)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sequence_id", "sequence", "sum_masked_log_likelihood", "mean_masked_log_probability", "position_log_probabilities"),
        )
        writer.writeheader()
        for identifier, sequence in read_fasta(input_fasta):
            position_scores = score_sequence(model, alphabet, sequence, device)
            writer.writerow({
                "sequence_id": identifier,
                "sequence": sequence,
                "sum_masked_log_likelihood": f"{sum(position_scores):.12g}",
                "mean_masked_log_probability": f"{sum(position_scores) / len(position_scores):.12g}",
                "position_log_probabilities": ";".join(f"{value:.12g}" for value in position_scores),
            })


def main() -> None:
    parser = argparse.ArgumentParser(description="ESM-2 masked pseudo-log-likelihood for peptide sequences")
    parser.add_argument("input_fasta", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--device", choices=("auto", "cpu", "cuda"), default="auto")
    args = parser.parse_args()
    run(args.input_fasta, args.output_csv, args.checkpoint, args.device)


if __name__ == "__main__":
    main()
