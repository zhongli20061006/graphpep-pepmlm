"""Pure identity and direction rules for the frozen CPU-only audit.

This module does not load models, infer experimental labels, normalize input
sequences, or calculate scores.  Near matches are HOLD records only.
"""

from __future__ import annotations

from collections.abc import Mapping


STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
NUMERICAL_CONTEXT = "FP32_eager_TF32off"
SCORE_DEFINITION = "masked_alt_minus_wt_logprob_v1"
P53_PARENT_13MER = "ETFSDLWKLLPEN"
PMI_PARENTS = frozenset(("TSFAEYWNLLSP", "TSFAEYWNLASP"))


def disagreement_type(dp: str, hp: str) -> str:
    """Separate opposing directions, abstention, ties, and method failures."""
    directional = {"gain", "loss"}
    if dp not in {"gain", "loss", "tie"}:
        return "method_failure"
    if hp not in {"gain", "loss", "tie", "abstain"}:
        return "method_failure"
    if dp in directional and hp in directional:
        return "same_direction" if dp == hp else "strict_opposite"
    if dp in directional:
        return "H_abstains"
    if hp in directional:
        return "d_tie_H_directional"
    return "both_nondirectional"


def _pair_sequences(pair: Mapping[str, object]) -> tuple[str, str] | None:
    """Read nonempty complete sequence fields without changing their content."""
    parent = pair.get("parent_sequence")
    variant = pair.get("variant_sequence")
    if not isinstance(parent, str) or not parent:
        return None
    if not isinstance(variant, str) or not variant:
        return None
    return parent, variant


def exact_pair(a: Mapping[str, object], b: Mapping[str, object]) -> bool:
    """Require both complete peptide endpoints and any supplied receptor."""
    a_sequences = _pair_sequences(a)
    b_sequences = _pair_sequences(b)
    if a_sequences is None or a_sequences != b_sequences:
        return False
    if "receptor_sequence" in a or "receptor_sequence" in b:
        a_receptor = a.get("receptor_sequence")
        b_receptor = b.get("receptor_sequence")
        return (
            isinstance(a_receptor, str)
            and bool(a_receptor)
            and isinstance(b_receptor, str)
            and a_receptor == b_receptor
        )
    return True


def variant_sequence(parent: str, site: int, alt: str) -> str:
    """Construct an identity string for one 1-based standard-AA replacement."""
    if not isinstance(parent, str) or not parent:
        raise ValueError("parent must be a nonempty standard-amino-acid string")
    if any(residue not in STANDARD_AMINO_ACIDS for residue in parent):
        raise ValueError("parent contains a nonstandard amino acid")
    if not isinstance(site, int) or isinstance(site, bool):
        raise ValueError("site must be a 1-based integer")
    if not 1 <= site <= len(parent):
        raise ValueError("site is outside the complete parent sequence")
    if not isinstance(alt, str) or len(alt) != 1:
        raise ValueError("alt must be one standard amino acid")
    if alt not in STANDARD_AMINO_ACIDS:
        raise ValueError("alt must be one standard amino acid")
    return parent[: site - 1] + alt + parent[site:]


def distribution_key(
    weight: str, receptor_sequence: str, parent_sequence: str, site: int
) -> tuple[str, str, str, int, str, str]:
    """Keep the complete model, sequence, site, and frozen score context."""
    return (
        weight,
        receptor_sequence,
        parent_sequence,
        site,
        NUMERICAL_CONTEXT,
        SCORE_DEFINITION,
    )


def _single_substitution(parent: str, variant: str) -> tuple[int, str, str] | None:
    """Derive a mutation from full endpoints, ignoring reported mutation names."""
    if len(parent) != len(variant):
        return None
    changes = [
        (site, wild_type, alternative)
        for site, (wild_type, alternative) in enumerate(zip(parent, variant), 1)
        if wild_type != alternative
    ]
    return changes[0] if len(changes) == 1 else None


def near_pair_reason(a: Mapping[str, object], b: Mapping[str, object]) -> str:
    """Identify only the specified HOLD relationships; never grant a label."""
    a_sequences = _pair_sequences(a)
    b_sequences = _pair_sequences(b)
    if a_sequences is None or b_sequences is None:
        return ""
    a_parent, a_variant = a_sequences
    b_parent, b_variant = b_sequences
    if a_sequences == b_sequences:
        return ""

    for long_pair, short_pair in (
        (a_sequences, b_sequences),
        (b_sequences, a_sequences),
    ):
        long_parent, long_variant = long_pair
        short_parent, short_variant = short_pair
        if (
            long_parent == P53_PARENT_13MER
            and len(short_parent) == 12
            and long_parent[:12] == short_parent
            and len(long_variant) == 13
            and len(short_variant) == 12
            and long_variant[:12] == short_variant
        ):
            return "HOLD_P53_13MER_VS_12MER"

    if {a_parent, b_parent} == PMI_PARENTS:
        a_mutation = _single_substitution(a_parent, a_variant)
        b_mutation = _single_substitution(b_parent, b_variant)
        if a_mutation is not None and a_mutation == b_mutation:
            return "HOLD_PMI_L10A_BACKGROUND_MISMATCH"

    if a_parent != b_parent and a_variant == b_variant:
        return "HOLD_SHARED_VARIANT_DIFFERENT_PARENT"
    return ""
