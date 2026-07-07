"""Effect-size and standard-error helpers for meta-analysis.

When a study does not report a standard error directly, we back it out from
the reported effect size and p-value using the normal approximation
se = |effect| / |z|, where z = Phi^-1(1 - p/2). This is a standard, documented
approximation (not a silent guess) and is only used when standard_error is
missing.
"""
from __future__ import annotations

from scipy.stats import norm

MIN_P = 1e-300


def approximate_se(effect_size: float, p_value: float) -> float | None:
    if effect_size is None or p_value is None:
        return None
    p = max(min(p_value, 1.0 - 1e-16), MIN_P)
    if p >= 1.0:
        return None
    z = norm.ppf(1 - p / 2)
    if z <= 0 or effect_size == 0:
        return None
    return abs(effect_size) / z


def effective_standard_error(effect_size: float, standard_error, p_value) -> float | None:
    if standard_error is not None and standard_error == standard_error and standard_error > 0:
        return float(standard_error)
    return approximate_se(effect_size, p_value)
