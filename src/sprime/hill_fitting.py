"""
Hill curve fitting module for sprime.

Implements a **linear-x four-parameter logistic (linear-x 4PL)** model: the independent
variable is concentration x on a linear scale, entering as (x / EC50) in the denominator.
This matches forms such as AAT Bioquest and generic ``(x/C)^B`` calculators.

Industry note: "Hill equation" and "4PL" are not standardized. Many tools (e.g. SigmaPlot,
GraphPad-style log-dose forms) use a **log-x 4PL** where the dose axis is log10(concentration),
which yields different Hill slope sign conventions for the same curve. Our slope values are
numerically consistent with the linear-x parameterization; compare log-x vs linear-x in
``docs/background/README_4PL_Dose_Response.md#linear-x-vs-log-x-4pl-hill-slope``.

**Naming:** The exponent *n* in ``(x/EC50)^n`` is exposed as ``steepness_coefficient`` (and
``initial_steepness_coefficient`` for guesses). That is the same role as the classical **Hill
coefficient** *n* in linear-x Hill formulations; we avoid the name ``hill_coefficient`` here
so it is not confused with log-x "Hill slope" outputs from other packages.

**Sign convention:** After fitting, parameters are always returned in **canonical form** using
Kendall concordance over all dose-response pairs to detect the biological curve direction. For
inhibitory curves (response decreases with dose), this guarantees ``zero_asymptote > inf_asymptote``
and therefore S' > 0. For disinhibitory curves, ``zero_asymptote < inf_asymptote`` and S' < 0. This
prevents the degenerate negative-*n* parameterization that ``scipy`` can converge to for steep
inhibitory curves, which would otherwise silently flip the sign of S'.

Adapted to work with sprime's domain entities.
"""

from typing import TYPE_CHECKING, List, Optional, Tuple

if TYPE_CHECKING:
    import numpy as np

# Import scipy (numpy comes as scipy dependency)
# Import at module level for type hints, but handle ImportError gracefully
try:
    import numpy as np
    from scipy.optimize import curve_fit
except ImportError:
    np = None
    curve_fit = None


def hill_equation(
    x, zero_asymptote: float, steepness_coefficient: float, ec50: float, inf_asymptote: float
):
    """
    Linear-x four-parameter logistic (linear-x 4PL) Hill form.

    Response vs concentration x (linear scale, same units as EC50):

        y = inf_asymptote + (zero_asymptote - inf_asymptote) / (1 + (x / C)^n)

    At concentration approaching zero, y approaches zero_asymptote; at saturating concentration,
    y approaches inf_asymptote. C = ec50, n = steepness_coefficient. Signed n encodes curve direction
    together with asymptote ordering; this is not interchangeable with log-x 4PL slope reporting.

    Args:
        x: Concentration values (linear scale)
        zero_asymptote: Response as concentration -> 0 (left side of dose axis)
        steepness_coefficient: Exponent n in (x/C)^n. Conceptually the same quantity often called
            the **Hill coefficient** in linear-x dose-response formulations (cooperativity exponent);
            we use ``steepness_coefficient`` here to avoid confusion with log-x tools that report a
            different "Hill slope" (see module docstring).
        ec50: Half-maximal concentration C
        inf_asymptote: Response at saturating concentration (right of dose axis)

    Returns:
        Response values
    """
    return inf_asymptote + (zero_asymptote - inf_asymptote) / (
        1 + (x / ec50) ** steepness_coefficient
    )


def _direction_is_down(responses) -> bool:
    """
    Detect inhibitory (decreasing) curve direction via Kendall concordance over all pairs.
    Ties default to inhibitory — conservative for drug screening where most active compounds
    reduce viability.
    """
    n = len(responses)
    down = sum(1 for i in range(n) for j in range(i + 1, n) if responses[j] < responses[i])
    up = sum(1 for i in range(n) for j in range(i + 1, n) if responses[j] > responses[i])
    return down >= up


def _normalize_direction(params, y_data):
    """
    Enforce canonical parameterization after fitting:
      inhibitory (data goes down) → zero_asymptote > inf_asymptote → S' positive
      disinhibitory  (data goes up)  → zero_asymptote < inf_asymptote → S' negative

    Swapping zero↔inf and flipping sign(n) describes the identical curve — same residuals,
    same r². This corrects the degenerate negative-n solution that scipy can converge to
    for steep inhibitory curves when seeded with n=-0.3.
    """
    from .sprime import HillCurveParams

    goes_down = _direction_is_down(list(y_data))
    zero_gt_inf = params.zero_asymptote > params.inf_asymptote

    if goes_down == zero_gt_inf:
        # Parameterization already matches the data direction:
        #   inhibitory (goes_down=True)  with zero > inf  -> S' positive (correct)
        #   disinhibitory (goes_down=False) with zero < inf  -> S' negative (correct)
        return params

    # Parameterization is in the degenerate form: asymptote ordering contradicts the data
    # direction.  Swap zero<->inf and flip sign(n) to recover the canonical form.
    # The two parameterizations describe the same curve -- r_squared is preserved.
    #   degenerate inhibitory: goes_down=True,  zero < inf  -> swap gives zero > inf, S' positive
    #   degenerate disinhibitory: goes_down=False, zero > inf  -> swap gives zero < inf, S' negative
    return HillCurveParams(
        ec50=params.ec50,
        zero_asymptote=params.inf_asymptote,
        inf_asymptote=params.zero_asymptote,
        steepness_coefficient=-params.steepness_coefficient,
        r_squared=params.r_squared,
    )


def fit_hill_curve(
    concentrations: List[float],
    responses: List[float],
    *,
    # Initial parameter guesses (all optional with defaults)
    initial_zero_asymptote: Optional[float] = None,
    initial_inf_asymptote: Optional[float] = None,
    initial_ec50: Optional[float] = None,
    initial_steepness_coefficient: Optional[float] = None,
    # Curve direction
    curve_direction: Optional[str] = None,  # "up", "down", or None for auto-detect
    # Optimization parameters
    maxfev: int = 3000000,
    # Zero concentration handling
    zero_replacement: float = 1e-24,
    # Parameter bounds (optional)
    bounds: Optional[Tuple[List[float], List[float]]] = None,
    # Additional scipy.optimize.curve_fit parameters
    **curve_fit_kwargs,
):
    """
    Fit four-parameter Hill equation to dose-response data.

    Fits a sigmoidal curve to concentration-response data and returns
    HillCurveParams with fitted parameters.

    Args:
        concentrations: List of concentration values
        responses: List of response values (must match length of concentrations)
        initial_zero_asymptote: Initial guess for zero asymptote (default: auto-estimated)
        initial_inf_asymptote: Initial guess for inf asymptote (default: auto-estimated)
        initial_ec50: Initial guess for EC50 (default: auto-estimated)
        initial_steepness_coefficient: Initial guess for steepness coefficient n (default: auto-estimated)
        curve_direction: Curve direction - "up" (increasing), "down" (decreasing),
                        or None for auto-detect (tries both, selects best r-squared)
        maxfev: Maximum function evaluations for optimization (default: 3,000,000)
        zero_replacement: Value to replace zero concentrations (default: 1e-24)
        bounds: Optional parameter bounds as (lower_bounds, upper_bounds) tuples
                Format: ([zero_asymptote_min, steepness_min, ec50_min, inf_asymptote_min],
                        [zero_asymptote_max, steepness_max, ec50_max, inf_asymptote_max])
        **curve_fit_kwargs: Additional arguments passed to scipy.optimize.curve_fit

    Returns:
        HillCurveParams: Fitted curve parameters. steepness_coefficient sign and thus S' direction
        are determined by Kendall concordance: inhibitory curves (positive S') or disinhibitory
        curves (negative S').

    Raises:
        ValueError: If inputs are invalid
        RuntimeError: If curve fitting fails
        ImportError: If numpy/scipy are not installed
    """
    if np is None or curve_fit is None:
        raise ImportError("Hill curve fitting requires scipy. " "Install with: pip install scipy")

    # Validate inputs
    if len(concentrations) != len(responses):
        raise ValueError("Concentrations and responses must have same length")

    if len(concentrations) < 4:
        raise ValueError("Need at least 4 data points to fit 4-parameter Hill equation")

    # Convert to numpy arrays (make copies to avoid modifying originals)
    concentrations = list(concentrations)
    responses = list(responses)

    # Sort data if needed (ascending concentrations)
    if concentrations[0] > concentrations[-1]:
        concentrations.reverse()
        responses.reverse()

    # Handle zero concentrations
    if concentrations[0] == 0:
        concentrations[0] = zero_replacement

    x_data = np.array(concentrations)
    y_data = np.array(responses)

    # Auto-detect or use specified curve direction
    if curve_direction is None:
        result = _fit_with_auto_direction(
            x_data,
            y_data,
            initial_zero_asymptote,
            initial_inf_asymptote,
            initial_ec50,
            initial_steepness_coefficient,
            maxfev,
            bounds,
            **curve_fit_kwargs,
        )
    else:
        result = _fit_single_direction(
            x_data,
            y_data,
            curve_direction,
            initial_zero_asymptote,
            initial_inf_asymptote,
            initial_ec50,
            initial_steepness_coefficient,
            maxfev,
            bounds,
            **curve_fit_kwargs,
        )

    # Enforce canonical sign convention regardless of which scipy basin was found.
    # Inhibitory (data goes down) -> zero > inf -> S' positive.
    # Activatory (data goes up)   -> zero < inf -> S' negative.
    return _normalize_direction(result, y_data)


def _fit_single_direction(
    x_data: "np.ndarray",
    y_data: "np.ndarray",
    curve_direction: str,
    initial_zero_asymptote: Optional[float],
    initial_inf_asymptote: Optional[float],
    initial_ec50: Optional[float],
    initial_steepness_coefficient: Optional[float],
    maxfev: int,
    bounds: Optional[Tuple[List[float], List[float]]],
    **curve_fit_kwargs,
):
    """Fit curve with specified direction."""
    # Import here to avoid circular import
    from .sprime import HillCurveParams

    # Get initial guesses (use defaults if not provided)
    if curve_direction == "up":
        # For curves that go up (increasing response)
        guess_zero_asymptote = (
            initial_zero_asymptote if initial_zero_asymptote is not None else 0.001
        )
        guess_steepness = (
            initial_steepness_coefficient if initial_steepness_coefficient is not None else 1.515
        )
        guess_ec50 = initial_ec50 if initial_ec50 is not None else 108.0
        guess_inf_asymptote = initial_inf_asymptote if initial_inf_asymptote is not None else 3.784
    else:  # "down"
        # For curves that go down (decreasing response)
        guess_zero_asymptote = (
            initial_zero_asymptote if initial_zero_asymptote is not None else 10.0
        )
        guess_steepness = (
            initial_steepness_coefficient if initial_steepness_coefficient is not None else -0.3
        )
        guess_ec50 = initial_ec50 if initial_ec50 is not None else 0.4
        guess_inf_asymptote = initial_inf_asymptote if initial_inf_asymptote is not None else 90.0

    # If any parameter not provided, try to estimate from data
    if initial_zero_asymptote is None:
        guess_zero_asymptote = (
            min(y_data) * 0.1 if guess_zero_asymptote == 0.001 else guess_zero_asymptote
        )
    if initial_inf_asymptote is None:
        guess_inf_asymptote = (
            max(y_data) * 1.1 if guess_inf_asymptote in (3.784, 90.0) else guess_inf_asymptote
        )
    if initial_ec50 is None:
        # Estimate EC50 as median concentration
        guess_ec50 = float(np.median(x_data))

    initial_guess = [guess_zero_asymptote, guess_steepness, guess_ec50, guess_inf_asymptote]

    # Fit the 4PL model
    try:
        fit_kwargs = {"p0": initial_guess, "maxfev": maxfev, **curve_fit_kwargs}
        if bounds is not None:
            fit_kwargs["bounds"] = bounds

        params, covariance = curve_fit(hill_equation, x_data, y_data, **fit_kwargs)
    except Exception as e:
        raise RuntimeError(f"Failed to fit Hill curve: {e}") from e

    # Extract fitted parameters
    zero_asymptote_fit, steepness_fit, ec50_fit, inf_asymptote_fit = params

    # Validate results
    if np.isnan(zero_asymptote_fit) or np.isnan(ec50_fit) or np.isnan(inf_asymptote_fit):
        raise RuntimeError(
            f"Fitting produced invalid parameters: zero_asymptote={zero_asymptote_fit}, "
            f"ec50={ec50_fit}, inf_asymptote={inf_asymptote_fit}"
        )

    # Calculate r-squared
    y_pred = hill_equation(x_data, *params)
    residuals = y_data - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return HillCurveParams(
        ec50=float(ec50_fit),
        zero_asymptote=float(zero_asymptote_fit),
        inf_asymptote=float(inf_asymptote_fit),
        steepness_coefficient=float(steepness_fit),
        r_squared=float(r_squared),
    )


def _fit_with_auto_direction(
    x_data: "np.ndarray",
    y_data: "np.ndarray",
    initial_zero_asymptote: Optional[float],
    initial_inf_asymptote: Optional[float],
    initial_ec50: Optional[float],
    initial_steepness_coefficient: Optional[float],
    maxfev: int,
    bounds: Optional[Tuple[List[float], List[float]]],
    **curve_fit_kwargs,
):
    """Try both curve directions and return best fit (highest r-squared)."""
    # Import here to avoid circular import
    best_params = None
    best_r2 = None

    for direction in ["up", "down"]:
        try:
            params = _fit_single_direction(
                x_data,
                y_data,
                direction,
                initial_zero_asymptote,
                initial_inf_asymptote,
                initial_ec50,
                initial_steepness_coefficient,
                maxfev,
                bounds,
                **curve_fit_kwargs,
            )

            if params.r_squared is not None:
                if best_r2 is None or params.r_squared > best_r2:
                    best_params = params
                    best_r2 = params.r_squared
        except (RuntimeError, ValueError):
            # Try next direction if this one fails
            continue

    if best_params is None:
        raise RuntimeError(
            "Failed to fit Hill curve in either direction. "
            "Check data quality and initial parameter guesses."
        )

    return best_params
