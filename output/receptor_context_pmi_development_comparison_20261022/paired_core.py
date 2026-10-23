"""Model-free, fixed-rule helpers for the finite PMI development comparison.

Sensitivity ranges are approximate delta-method scenarios, not source confidence
intervals or guaranteed conservative confidence intervals.  Summary counts are
descriptive development results; this module does not calculate formal BA or CI.
"""
from __future__ import annotations

from decimal import Decimal, InvalidOperation
from html.parser import HTMLParser
import math
import re
from typing import Any, Iterable, Mapping


EPSILON = 1e-4
T_DF2 = 4.302652729911275
N_LOWER_BOUND_ASSUMED = 3
SCORING_MARKER = "single_mask_alt_minus_wt_logprob_v1"

_MANTISSA = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)"
_KD_PATTERN = re.compile(
    rf"^({_MANTISSA})(?:[eE]([+-]?\d+))?"
    rf"\s*±\s*({_MANTISSA})(?:[eE]([+-]?\d+))?$"
)


def parse_kd(text: str) -> tuple[float, float] | None:
    """Convert a molar mean ± SD to nM, including a shared final exponent.

    In ``3.2±1.1E-09``, the exponent applies to both printed mantissas.
    Decimal scaling occurs before conversion to float, avoiding scale rounding.
    """
    value = " ".join(str(text).replace("−", "-").split())
    missing = re.sub(r"[.\s]", "", value).lower()
    if missing in {"nd", "nda"}:
        return None
    match = _KD_PATTERN.fullmatch(value)
    if match is None:
        raise ValueError(f"Unsupported Kd mean/SD notation: {text!r}")
    mean_text, mean_exp, sd_text, sd_exp = match.groups()
    # A lone trailing exponent is shared; explicit exponents remain separate.
    if mean_exp is None and sd_exp is not None:
        mean_exp = sd_exp
    try:
        mean_nm = Decimal(mean_text).scaleb(int(mean_exp or 0) + 9)
        sd_nm = Decimal(sd_text).scaleb(int(sd_exp or 0) + 9)
        mean, sd = float(mean_nm), float(sd_nm)
    except (InvalidOperation, OverflowError, ValueError) as exc:
        raise ValueError(f"Invalid Kd mean/SD: {text!r}") from exc
    _validate_measurement(mean, sd)
    return mean, sd


class _TableParser(HTMLParser):
    """Read td data rows; rowspan intentionally does not shift source columns."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.tables: list[list[list[str]]] = []
        self.table: list[list[str]] | None = None
        self.row: list[str] | None = None
        self.cell_parts: list[str] | None = None
        self.colspan = 1

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "table":
            self.table = []
        elif tag == "tr" and self.table is not None:
            self.row = []
        elif tag == "td" and self.row is not None:
            self.cell_parts = []
            try:
                self.colspan = max(1, int(dict(attrs).get("colspan") or 1))
            except ValueError:
                self.colspan = 1

    def handle_data(self, data: str) -> None:
        if self.cell_parts is not None:
            self.cell_parts.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag == "td" and self.cell_parts is not None and self.row is not None:
            cell = " ".join("".join(self.cell_parts).split())
            self.row.extend([cell] + [""] * (self.colspan - 1))
            self.cell_parts = None
        elif tag == "tr" and self.row is not None:
            if self.row and self.table is not None:
                self.table.append(self.row)
            self.row = None
        elif tag == "table" and self.table is not None:
            self.tables.append(self.table)
            self.table = None


def parse_tables(html: str) -> list[list[list[str]]]:
    """Read each non-nested HTML table, ignoring th and expanding td colspan."""
    parser = _TableParser()
    parser.feed(html)
    parser.close()
    return parser.tables


def _validate_measurement(mean: float, sd: float = 0.0) -> None:
    if not math.isfinite(mean) or mean <= 0:
        raise ValueError("Kd mean must be finite and positive")
    if not math.isfinite(sd) or sd < 0:
        raise ValueError("Kd SD must be finite and nonnegative")


def point_direction(variant_mean: float, parent_mean: float) -> str:
    """Apply the fixed twofold rule to the unrounded mean ratio."""
    _validate_measurement(variant_mean)
    _validate_measurement(parent_mean)
    ratio = variant_mean / parent_mean
    if ratio <= 0.5:
        return "gain"
    if ratio >= 2.0:
        return "loss"
    return "ambiguous"


def sensitivity(
    variant_mean: float, variant_sd: float, parent_mean: float, parent_sd: float
) -> dict[str, Any]:
    """Approximate log-ratio sensitivity at assumed n=3 and t(df=2).

    Independent errors combine variances; the worst-correlation scenario adds
    the two relative SDs.  Exact replicate counts and raw repeats are unknown.
    """
    _validate_measurement(variant_mean, variant_sd)
    _validate_measurement(parent_mean, parent_sd)
    ratio = variant_mean / parent_mean
    log_ratio = math.log(variant_mean) - math.log(parent_mean)
    relative_variant = variant_sd / variant_mean
    relative_parent = parent_sd / parent_mean
    root_n = math.sqrt(N_LOWER_BOUND_ASSUMED)
    se_independent = math.hypot(relative_variant, relative_parent) / root_n
    se_worst = (relative_variant + relative_parent) / root_n
    return {
        "ratio": ratio,
        "independent_lower": math.exp(log_ratio - T_DF2 * se_independent),
        "independent_upper": math.exp(log_ratio + T_DF2 * se_independent),
        "worst_correlation_lower": math.exp(log_ratio - T_DF2 * se_worst),
        "worst_correlation_upper": math.exp(log_ratio + T_DF2 * se_worst),
        "n_lower_bound_assumed": N_LOWER_BOUND_ASSUMED,
        "t_df2": T_DF2,
        "se_log_ratio_independent": se_independent,
        "se_log_ratio_worst_correlation": se_worst,
        "method": "approximate_log_ratio_delta_method_sensitivity",
        "is_source_confidence_interval": False,
        "is_guaranteed_conservative_confidence_interval": False,
    }


def predict_d(value: float) -> str:
    """Predict a direction with the frozen 1e-4-nat tie tolerance."""
    if not math.isfinite(value):
        return "failure"
    if value > EPSILON:
        return "gain"
    if value < -EPSILON:
        return "loss"
    return "tie"


def predict_h(control_1: float, control_2: float) -> str:
    """Require both fixed control corrections to predict the same direction."""
    if not math.isfinite(control_1) or not math.isfinite(control_2):
        return "failure"
    if control_1 > EPSILON and control_2 > EPSILON:
        return "gain"
    if control_1 < -EPSILON and control_2 < -EPSILON:
        return "loss"
    return "abstain"


_COUNT_KEYS = (
    "n", "d_correct", "h_correct", "both_correct", "h_only_correct",
    "d_only_correct", "neither_correct", "d_tie_or_failure",
    "h_abstain_or_failure",
)


def paired_summary(rows: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    """Count paired outcomes over every supplied gain/loss development pair.

    Abstentions, ties, missing predictions and failures remain in all applicable
    denominators.  Unknown/ambiguous reference labels must be audited separately.
    """
    result: dict[str, Any] = dict.fromkeys(_COUNT_KEYS, 0)
    by_direction = {
        direction: dict.fromkeys(_COUNT_KEYS, 0) for direction in ("gain", "loss")
    }
    for row in rows:
        direction = row["point_direction"]
        if direction not in by_direction:
            raise ValueError("paired_summary requires gain/loss reference labels")
        d_prediction = row.get("d_prediction", "failure")
        h_prediction = row.get("h_prediction", "failure")
        d_correct = d_prediction == direction
        h_correct = h_prediction == direction
        increments = {
            "n": 1,
            "d_correct": int(d_correct),
            "h_correct": int(h_correct),
            "both_correct": int(d_correct and h_correct),
            "h_only_correct": int(h_correct and not d_correct),
            "d_only_correct": int(d_correct and not h_correct),
            "neither_correct": int(not d_correct and not h_correct),
            "d_tie_or_failure": int(d_prediction not in ("gain", "loss")),
            "h_abstain_or_failure": int(h_prediction not in ("gain", "loss")),
        }
        for key, amount in increments.items():
            result[key] += amount
            by_direction[direction][key] += amount
    result["by_direction"] = by_direction
    return result


def cache_key(
    checkpoint: str, receptor: str, peptide: str, site: int, precision: str
) -> tuple[str, str, str, int, str, str]:
    """Key a base distribution by its full exact context and fixed scoring rule."""
    return checkpoint, receptor, peptide, site, precision, SCORING_MARKER
