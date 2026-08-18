"""
Activity calculation from identified nuclides (Phase 2.1d).

The central formula is:

    A = S_net / (ε(E) · I_γ · t_live · p)

where:
    S_net   — net photopeak area on background-subtracted spectrum [counts]
    ε(E)    — photopeak efficiency at line energy [dimensionless]
    I_γ     — emission probability per decay [decimal, e.g. 0.851 for Cs-137]
    t_live  — live time [s]
    p       — correction factor (coincidence summing etc.); default 1.0

For nuclides with multiple matched lines, the per-line activities are
combined by weighted averaging with weights w_i = 1/σ²(A_i) (Gilmore
§5.7.2). The intra-nuclide χ²/dof reports how well the individual
lines agree — a poor χ² indicates either (a) cascade-summing depletion
on some lines, (b) an interfering nuclide contributing to one of the
peaks, or (c) the line was misassigned.

Reference: Gilmore G., Joss D. "Practical Gamma-ray Spectrometry"
3rd Ed., Wiley 2024, §5.7; Lsrm Algorithmic Foundations 2022, §8.4.
"""

from gamma.activity.compute import (
    BackgroundNotSubtractedError,
    LineActivity,
    ActivityResult,
    CASCADE_SUMMING_NUCLIDES,
    DEFAULT_TCS_METHOD_SCALE,
    compute_activity,
    compute_activities_for_all,
)

__all__ = [
    "BackgroundNotSubtractedError",
    "LineActivity",
    "ActivityResult",
    "CASCADE_SUMMING_NUCLIDES",
    "DEFAULT_TCS_METHOD_SCALE",
    "compute_activity",
    "compute_activities_for_all",
]
