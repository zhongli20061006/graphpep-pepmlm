"""Same-assay, same-construct MDM2 references from Li et al. 2010 Tables 1/3."""
from __future__ import annotations

from dataclasses import dataclass
import math
from statistics import median

from scipy.stats import spearmanr

SOURCE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC2856455/"
SOURCE_DOI = "10.1016/j.jmb.2010.03.005"


@dataclass(frozen=True)
class ReferencePeptide:
    peptide_id: str
    family: str
    sequence: str
    kd_nM: float
    sd_nM: float
    parent_id: str
    position: int
    mutation: str
    table: int


REFERENCES = (
    ReferencePeptide("PMI_WT", "PMI", "TSFAEYWNLLSP", 3.2, 1.1, "PMI_WT", 0, "WT", 1),
    ReferencePeptide("PMI_N8A", "PMI", "TSFAEYWALLSP", 0.49, 0.21, "PMI_WT", 8, "N8A", 1),
    ReferencePeptide("PMI_Y6A", "PMI", "TSFAEAWNLLSP", 610, 70, "PMI_WT", 6, "Y6A", 1),
    ReferencePeptide("PMI_F3A", "PMI", "TSAAEYWNLLSP", 38000, 2000, "PMI_WT", 3, "F3A", 1),
    ReferencePeptide("PMI_W7A", "PMI", "TSFAEYANLLSP", 160000, 30000, "PMI_WT", 7, "W7A", 1),
    ReferencePeptide("P53_WT", "p53", "ETFSDLWKLLPE", 440, 40, "P53_WT", 0, "WT", 3),
    ReferencePeptide("P53_P27A", "p53", "ETFSDLWKLLAE", 51, 3, "P53_WT", 11, "P27A", 3),
    ReferencePeptide("P53_K24A", "p53", "ETFSDLWALLPE", 230, 20, "P53_WT", 8, "K24A", 3),
    ReferencePeptide("P53_L22A", "p53", "ETFSDAWKLLPE", 5000, 400, "P53_WT", 6, "L22A", 3),
    ReferencePeptide("P53_L26A", "p53", "ETFSDLWKLAPE", 27000, 1000, "P53_WT", 10, "L26A", 3),
)


def affinity_rank_correlation(kd_nM: list[float], scores: list[float]) -> float:
    """Positive rho means lower GraphPep scores align with lower experimental Kd."""
    if len(kd_nM) != len(scores) or len(scores) < 3:
        raise ValueError("Need at least three paired peptide observations")
    if any(not math.isfinite(value) or value <= 0 for value in kd_nM):
        raise ValueError("Kd must be finite and positive")
    if any(not math.isfinite(value) for value in scores):
        raise ValueError("Scores must be finite")
    if len(set(scores)) < 2:
        return math.nan
    return float(spearmanr([math.log10(value) for value in kd_nM], scores).statistic)


def paired_median_gain(parent: list[float], child: list[float]) -> float:
    """Use paired differences, preserving AFM seed/model identity."""
    if len(parent) != len(child) or len(parent) < 2:
        raise ValueError("Need at least two matched models")
    if any(not math.isfinite(value) for value in parent + child):
        raise ValueError("Scores must be finite")
    return float(median(left - right for left, right in zip(parent, child)))
