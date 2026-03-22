"""
Response preprocessing for raw qHTS curves before linear-x 4PL fitting.

Implements the **test/control ratio**, **max normalization** (scale so the largest value is **1**,
then x100), and the conventional **x100** scale so that intermediate values align with the committed
test fixture ``tests/fixtures/SPrime_variation_reference.csv`` (columns *test/control ratio*,
*normalise to 1*, *x 100 scale*, and *Nonnormalized, x100 scale*).

For typical **inhibition-style** curves (response falls with dose), the largest control-relative
ratio is usually at the **lowest concentration**; then :func:`normalize_to_max_value` matches the
spreadsheet *normalise to 1* step (first point -> 1). If a later point is numerically larger, max
normalization scales to that peak instead of forcing index 0 to 1.

**DMSO** in docstrings means the **vehicle control** readout (typical solvent: DMSO). The
recommended CSV column for the per-profile control value at import time is
**``Control_Response``** (see project documentation).

This module is the **library** source of truth for the x100 convention; callers should use
these functions rather than re-implementing factors ad hoc. For **strict** CSV import
(``skip_control_response_normalization=False``), :meth:`~sprime.sprime.RawDataset.to_screening_dataset`
applies :func:`pipeline_asymptote_normalized` automatically before fitting raw curves.
Use these helpers directly when building in-memory arrays or custom pipelines (e.g. response-scale only).

The variation reference fixture is maintained so that *Nonnormalized, x100 scale* equals ``100 x`` the
second *test/control ratio* column row-by-row, and raw percent nucleus values are consistent
with those ratios given the DMSO control row (``35.3`` in that file). The library
definition ``pipeline_response_scale`` = ``100 x (raw / Control_Response)`` matches those
columns within normal floating-point tolerance when the same inputs are used.
"""

from __future__ import annotations

from typing import List, Sequence

#: Fixed multiplier (**100**) for both major raw pipelines (reference CSV *x 100 scale*).
#: Treat as **immutable in practice** for production - ``scale_factor`` / ``factor`` arguments
#: that default to this value should **stay at 100** unless there is a rare, documented exception.
S_PRIME_RESPONSE_SCALE_FACTOR: float = 100.0


def ratios_to_control(responses: Sequence[float], control_response: float) -> List[float]:
    """
    First preprocessing step: **test/control ratio** - divide each raw response by the
    vehicle (**DMSO**) control for that experimental context (e.g. same plate).

    Args:
        responses: Raw readouts at each concentration (same order as concentrations).
        control_response: Single control (**DMSO**) response; must be non-zero.

    Returns:
        List of ratios ``r / control_response``.

    Raises:
        ValueError: If ``control_response`` is zero.
    """
    if control_response == 0:
        raise ValueError(
            "control_response must be non-zero for control-relative ratios (DMSO/vehicle readout)."
        )
    cr = float(control_response)
    return [float(r) / cr for r in responses]


def normalize_to_max_value(values: Sequence[float]) -> List[float]:
    """
    Scale a series so the **maximum** value becomes **1.0**; all other values are <= 1.

    After **test/control** ratios, this is the asymptote-normalized path's "normalise to 1"
    step before x100: the **largest** ratio (often the low-dose / vehicle-region point on
    descending curves) maps to 1; remaining points are proportional and lower.

    Args:
        values: Non-empty sequence (e.g. control-relative ratios).

    Returns:
        List where ``max(result) == 1.0`` (within float tolerance) and each entry is
        ``values[i] / max(values)``.

    Raises:
        ValueError: If ``values`` is empty, or if ``max(values) == 0``.
    """
    if not values:
        raise ValueError("normalize_to_max_value requires a non-empty sequence")
    v = [float(x) for x in values]
    m = max(v)
    if m == 0:
        raise ValueError("normalize_to_max_value: maximum value is zero (cannot scale to unity).")
    return [x / m for x in v]


def scale_responses(
    values: Sequence[float],
    factor: float = S_PRIME_RESPONSE_SCALE_FACTOR,
) -> List[float]:
    """
    Apply the conventional **x100** scale.

    ``factor`` defaults to **100** and should **almost never** be changed for production
    workflows (matches the variation reference and S' conventions). Non-default values are for
    rare cases such as unit tests or hand-built pipelines.

    The factor is a **reporting convention** only; it does not add biological meaning by
    itself. Both *asymptote-normalized* and *response-scale* pipelines use this after the
    prior steps defined in the reference spreadsheet.
    """
    f = float(factor)
    return [float(x) * f for x in values]


def pipeline_response_scale(
    raw_responses: Sequence[float],
    control_response: float,
    *,
    scale_factor: float = S_PRIME_RESPONSE_SCALE_FACTOR,
) -> List[float]:
    """
    **Response-scale** path (no max normalization): test/control ratio, then multiply by ``scale_factor``.

    Matches the reference CSV path: *test/control ratio* -> *Nonnormalized, x100 scale*
    when ``scale_factor`` is ``100`` (the default - keep it at **100** unless you have a rare,
    documented reason to override).
    """
    ratios = ratios_to_control(raw_responses, control_response)
    return scale_responses(ratios, scale_factor)


def pipeline_asymptote_normalized(
    raw_responses: Sequence[float],
    control_response: float,
    *,
    scale_factor: float = S_PRIME_RESPONSE_SCALE_FACTOR,
) -> List[float]:
    """
    **Asymptote-normalized** path: test/control ratio -> :func:`normalize_to_max_value` -> multiply by ``scale_factor``.

    Matches the reference CSV columns *test/control ratio*, *normalise to 1*, *x 100 scale* when
    the max ratio coincides with the sheet's intended reference (typically the low-dose peak).

    ``scale_factor`` defaults to **100**; it should **rarely** differ from that for real data.
    """
    ratios = ratios_to_control(raw_responses, control_response)
    unity = normalize_to_max_value(ratios)
    return scale_responses(unity, scale_factor)
