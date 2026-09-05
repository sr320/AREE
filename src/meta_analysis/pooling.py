"""Random-effects (DerSimonian-Laird) meta-analysis, implemented directly for
auditability (see docs/design.md Assumptions — no compiled R meta package
dependency in the MVP), plus Benjamini-Hochberg false-discovery-rate control
for the family of pooled tests a run produces.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm


@dataclass
class PooledResult:
    k: int
    pooled_effect: float
    pooled_se: float
    ci_lower: float
    ci_upper: float
    z: float
    p_value: float
    q_statistic: float
    i_squared: float
    tau_squared: float


def two_sided_p_from_z(z: float) -> float:
    """Two-sided normal p-value that stays strictly positive for large |z|.

    `2 * (1 - norm.cdf(|z|))` underflows to exactly 0.0 once |z| exceeds ~8.3,
    which on a genome-wide pool happened for hundreds of genes and printed as
    "p = 0" on their evidence cards. The survival function keeps precision in
    the tail, so the ordering among very small p-values is preserved and the
    BH adjustment downstream sees real numbers.
    """
    return float(2.0 * norm.sf(abs(z)))


def dersimonian_laird(effect_sizes, standard_errors) -> PooledResult:
    """Pool a list of study-level effect sizes and standard errors.

    Falls back gracefully to a fixed-effect (single-study) result when k == 1.
    """
    yi = np.asarray(effect_sizes, dtype=float)
    sei = np.asarray(standard_errors, dtype=float)
    k = len(yi)
    if k == 0:
        raise ValueError("Cannot pool zero studies")

    wi = 1.0 / (sei ** 2)
    fixed_effect = float(np.sum(wi * yi) / np.sum(wi))

    if k == 1:
        pooled_se = float(sei[0])
        q_statistic = 0.0
        tau_squared = 0.0
        i_squared = 0.0
        pooled_effect = fixed_effect
    else:
        q_statistic = float(np.sum(wi * (yi - fixed_effect) ** 2))
        df = k - 1
        c = float(np.sum(wi) - np.sum(wi ** 2) / np.sum(wi))
        tau_squared = max(0.0, (q_statistic - df) / c) if c > 0 else 0.0
        wi_star = 1.0 / (sei ** 2 + tau_squared)
        pooled_effect = float(np.sum(wi_star * yi) / np.sum(wi_star))
        pooled_se = float(np.sqrt(1.0 / np.sum(wi_star)))
        i_squared = max(0.0, (q_statistic - df) / q_statistic) * 100 if q_statistic > 0 else 0.0

    z = pooled_effect / pooled_se if pooled_se > 0 else 0.0
    p_value = two_sided_p_from_z(z)
    ci_lower = pooled_effect - 1.96 * pooled_se
    ci_upper = pooled_effect + 1.96 * pooled_se

    return PooledResult(
        k=k,
        pooled_effect=pooled_effect,
        pooled_se=pooled_se,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        z=z,
        p_value=p_value,
        q_statistic=q_statistic,
        i_squared=i_squared,
        tau_squared=tau_squared,
    )


def benjamini_hochberg(p_values) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (FDR q-values) for one test family.

    Returns an array aligned with the input. NaN inputs stay NaN and do not
    count toward the family size. Values are monotone in the input ordering
    and clipped to [0, 1], following the standard step-up procedure.
    """
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = ~np.isnan(p)
    n = int(valid.sum())
    if n == 0:
        return out
    pv = p[valid]
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    restored = np.empty(n)
    restored[order] = adjusted
    out[valid] = restored
    return out
