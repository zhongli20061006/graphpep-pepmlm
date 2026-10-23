"""Reuse a receptor-only alignment across peptide sequences and lengths."""
from __future__ import annotations

import hashlib

from feedback_pilot import aligned, parse_a3m


def make_shared_msa(text: str, receptor: str, source_peptide: str, peptide: str) -> tuple[str, dict]:
    """Replace the peptide query and padding without changing receptor alignment."""
    nr, ns, np = len(receptor), len(source_peptide), len(peptide)
    records = parse_a3m(text)
    if not ns or not np or not records or aligned(records[0][1]) != receptor + source_peptide:
        raise ValueError("Cached MSA query does not match the source construct")
    receptor_rows = []
    for header, sequence in records[1:]:
        clean = aligned(sequence)
        if len(clean) != nr + ns:
            raise ValueError("Unexpected source A3M aligned length")
        if any(c != "-" for c in clean[:nr]) and set(clean[nr:]) == {"-"}:
            if not sequence.endswith("-" * ns):
                raise ValueError("Unexpected insertion in the peptide padding")
            receptor_rows.append((header, sequence[:-ns]))
    if not receptor_rows:
        raise ValueError("No receptor-only MSA rows")
    invariant_body = "".join(f">{header}\n{sequence}\n" for header, sequence in receptor_rows)
    body = "".join(f">{header}\n{sequence}{'-' * np}\n" for header, sequence in receptor_rows)
    result = f"#{nr},{np}\t1,1\n>101\t102\n{receptor}{peptide}\n" + body + f">102\n{'-' * nr}{peptide}\n"
    stats = {"source_records": len(records), "receptor_only_records": len(receptor_rows),
             "discarded_non_receptor_records": len(records) - 1 - len(receptor_rows),
             "source_peptide_length": ns, "peptide_length": np,
             "receptor_alignment_sha256": hashlib.sha256(invariant_body.encode()).hexdigest(),
             "common_receptor_rows_sha256": hashlib.sha256(body.encode()).hexdigest()}
    return result, stats
