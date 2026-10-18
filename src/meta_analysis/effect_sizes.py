"""Effect-size and standard-error helpers for meta-analysis.

When a study does not report a standard error directly, we back it out from
the reported effect size and p-value using the normal approximation
se = |effect| / |z|, where z = Phi^-1(1 - p/2). This is a standard, documented
approximation (not a silent guess) and is only used when standard_error is
missing.
"""
from __future__ import annotations

import math

from scipy.stats import norm

MIN_P = 1e-300


def approximate_se(effect_size: float, p_value: float) -> float | None:
    if (effect_size is None or p_value is None
            or not math.isfinite(effect_size) or not math.isfinite(p_value)
            or not 0 <= p_value < 1 or effect_size == 0):
        return None
    p = max(p_value, MIN_P)
    # Subtracting p/2 from 1 rounds to 1 for very small p, giving z=inf
    # and a zero SE. The inverse survival function preserves the tail.
    z = norm.isf(p / 2)
    return abs(effect_size) / z


def effective_standard_error(effect_size: float, standard_error, p_value) -> float | None:
    if standard_error is not None and math.isfinite(standard_error) and standard_error > 0:
        return float(standard_error)
    return approximate_se(effect_size, p_value)
