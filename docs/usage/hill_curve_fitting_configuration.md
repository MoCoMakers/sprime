# Hill Curve Fitting Configuration Guide

This document describes all configurable parameters for Hill curve fitting in sprime. All parameters have sensible defaults but can be overridden to customize the fitting process.

## Overview

The Hill curve fitting uses a four-parameter logistic (4PL) model to fit dose-response data. All assumptions and defaults can be modified to suit your specific data characteristics.

## Function Signature

```python
hill_fitting.fit_hill_curve(
    concentrations: List[float],
    responses: List[float],
    *,
    initial_lower: Optional[float] = None,
    initial_upper: Optional[float] = None,
    initial_ec50: Optional[float] = None,
    initial_hill_coefficient: Optional[float] = None,
    curve_direction: Optional[str] = None,
    maxfev: int = 3000000,
    zero_replacement: float = 1e-24,
    bounds: Optional[Tuple[List[float], List[float]]] = None,
    **curve_fit_kwargs
) -> HillCurveParams
```

## Configurable Parameters

### Initial Parameter Guesses

These parameters provide starting values for the optimization algorithm. If not specified, defaults are used based on curve direction and data characteristics.

#### `initial_lower` (Optional[float], default: None)

Initial guess for the lower asymptote (minimum response value).

**Default behavior:**
- For "up" curves: `0.001` (or `min(responses) * 0.1` if auto-estimated)
- For "down" curves: `10.0` (or estimated from data)

**When to modify:**
- If your data has a known baseline response value
- If default guesses lead to poor convergence
- If responses are in a specific range (e.g., 0-100 for percentages)

**Example:**
```python
# If you know baseline is around 5%
hill_params = fit_hill_curve(
    concentrations, responses,
    initial_lower=5.0
)
```

#### `initial_upper` (Optional[float], default: None)

Initial guess for the upper asymptote (maximum response value).

**Default behavior:**
- For "up" curves: `3.784` (or `max(responses) * 1.1` if auto-estimated)
- For "down" curves: `90.0` (or estimated from data)

**When to modify:**
- If you know the maximum achievable response
- If default guesses are far from your data range
- For percentage data, might want to set to 100.0

**Example:**
```python
# For percentage response data
hill_params = fit_hill_curve(
    concentrations, responses,
    initial_upper=100.0
)
```

#### `initial_ec50` (Optional[float], default: None)

Initial guess for EC50 (half-maximal concentration).

**Default behavior:**
- For "up" curves: `108.0`
- For "down" curves: `0.4`
- If None: Uses `median(concentrations)` as estimate

**When to modify:**
- If you have prior knowledge of approximate EC50
- If default values are far from your concentration range
- Can significantly improve convergence speed

**Example:**
```python
# If you expect EC50 around 1 microM
hill_params = fit_hill_curve(
    concentrations, responses,
    initial_ec50=1.0
)
```

#### `initial_hill_coefficient` (Optional[float], default: None)

Initial guess for Hill coefficient (slope/steepness).

**Default behavior:**
- For "up" curves: `1.515`
- For "down" curves: `-0.3`

**When to modify:**
- If you know the expected steepness of your curve
- Typical values: 0.5 to 3.0 for most biological assays
- Negative values for decreasing curves

**Example:**
```python
# For a steep curve
hill_params = fit_hill_curve(
    concentrations, responses,
    initial_hill_coefficient=2.5
)
```

### Curve Direction

#### `curve_direction` (Optional[str], default: None)

Specifies whether the curve increases or decreases with concentration.

**Options:**
- `None`: Auto-detect (tries both "up" and "down", selects best R²)
- `"up"`: Increasing response with concentration (e.g., activation)
- `"down"`: Decreasing response with concentration (e.g., inhibition)

**Default behavior:**
- Auto-detects by trying both directions and selecting the fit with highest R²

**When to modify:**
- If you know the direction and want faster fitting (skip bidirectional)
- If auto-detection selects wrong direction
- For batch processing where direction is consistent

**Example:**
```python
# Force increasing curve
hill_params = fit_hill_curve(
    concentrations, responses,
    curve_direction="up"
)
```

### Optimization Parameters

#### `maxfev` (int, default: 3000000)

Maximum number of function evaluations during optimization.

**Default:** 3,000,000 (very high to ensure convergence)

**When to modify:**
- Reduce for faster fitting (e.g., `maxfev=10000`)
- Increase if convergence fails (though 3M is already very high)
- Typical values: 1000-10000 for most datasets

**Example:**
```python
# Faster fitting for simple curves
hill_params = fit_hill_curve(
    concentrations, responses,
    maxfev=10000
)
```

### Zero Concentration Handling

#### `zero_replacement` (float, default: 1e-24)

Value used to replace zero concentrations to avoid numerical issues.

**Default:** `1e-24` (extremely small)

**When to modify:**
- If you want a more reasonable small value (e.g., `1e-10`)
- If zero replacement affects your results
- Consider using bounds instead if zeros are problematic

**Example:**
```python
# Use a more reasonable small value
hill_params = fit_hill_curve(
    concentrations, responses,
    zero_replacement=1e-10
)
```

### Parameter Bounds

#### `bounds` (Optional[Tuple[List[float], List[float]]], default: None)

Constraints on parameter values during optimization.

**Format:**
```python
bounds = (
    [lower_min, hill_min, ec50_min, upper_min],  # Lower bounds
    [lower_max, hill_max, ec50_max, upper_max]   # Upper bounds
)
```

**Default:** None (unconstrained optimization)

**When to use:**
- To enforce biological constraints (e.g., EC50 must be positive)
- To prevent unrealistic parameter values
- To improve convergence when parameters are known to be in a range

**Example:**
```python
# Constrain EC50 to be between 0.1 and 100 microM
# Constrain response to be between 0 and 100%
hill_params = fit_hill_curve(
    concentrations, responses,
    bounds=(
        [0, -5, 0.1, 0],      # Lower bounds: [lower, hill, ec50, upper]
        [100, 5, 100, 100]     # Upper bounds
    )
)
```

### Additional scipy.optimize.curve_fit Parameters

Any additional keyword arguments are passed directly to `scipy.optimize.curve_fit`.

**Common options:**
- `method`: Optimization method (default: 'lm' for Levenberg-Marquardt)
- `absolute_sigma`: Whether to use absolute sigma values
- `check_finite`: Whether to check for finite values

**Example:**
```python
# Use different optimization method
hill_params = fit_hill_curve(
    concentrations, responses,
    method='trf'  # Trust Region Reflective algorithm
)
```

## Usage Examples

### Basic Usage (All Defaults)

```python
from sprime import hill_fitting

hill_params = hill_fitting.fit_hill_curve(
    concentrations=[0, 1, 10, 100, 1000],
    responses=[5, 10, 50, 90, 95]
)
```

### Custom Initial Guesses

```python
hill_params = hill_fitting.fit_hill_curve(
    concentrations=concentrations,
    responses=responses,
    initial_lower=0.0,      # Known baseline
    initial_upper=100.0,    # Percentage data
    initial_ec50=10.0,      # Approximate EC50
    initial_hill_coefficient=1.5
)
```

### Constrained Fitting

```python
hill_params = hill_fitting.fit_hill_curve(
    concentrations=concentrations,
    responses=responses,
    bounds=(
        [0, 0.5, 0.01, 0],      # Minimum values
        [100, 5.0, 1000, 100]   # Maximum values
    )
)
```

### Fast Fitting for Batch Processing

```python
hill_params = hill_fitting.fit_hill_curve(
    concentrations=concentrations,
    responses=responses,
    curve_direction="up",   # Skip auto-detection
    maxfev=5000            # Faster convergence
)
```

## Integration with DoseResponseProfile

When using `DoseResponseProfile.fit_hill_curve()`, all parameters are passed through:

```python
from sprime import DoseResponseProfile, Compound, CellLine, Assay

profile = DoseResponseProfile(
    compound=compound,
    cell_line=cell_line,
    assay=assay,
    concentrations=[...],
    responses=[...]
)

# Fit with custom parameters
hill_params = profile.fit_hill_curve(
    initial_ec50=1.0,
    maxfev=10000,
    curve_direction="up"
)
```

## Troubleshooting

### Convergence Failures

If fitting fails to converge:
1. **Check data quality**: Ensure sufficient data points (≥4) spanning full response range
2. **Adjust initial guesses**: Provide better starting values based on your data
3. **Increase maxfev**: Allow more iterations (though default is already very high)
4. **Add bounds**: Constrain parameters to reasonable ranges
5. **Try different curve direction**: Explicitly specify "up" or "down"

### Poor Fit Quality (Low R²)

If R² is low:
1. **Check data**: Ensure data follows sigmoidal pattern
2. **Remove outliers**: Outliers can significantly affect fit
3. **Adjust initial guesses**: Better starting values can improve convergence
4. **Try bounds**: Constraining parameters may help
5. **Check data range**: Ensure concentrations span EC50 region

### Invalid Parameters (NaN values)

If parameters are NaN:
1. **Check for zero/invalid concentrations**: Zero handling may need adjustment
2. **Verify data types**: Ensure all values are numeric
3. **Check data range**: Very large or very small values can cause issues
4. **Try different initial guesses**: Current guesses may be incompatible with data

## Default Values Summary

| Parameter | "Up" Curve Default | "Down" Curve Default | Auto-Estimated |
|-----------|-------------------|---------------------|----------------|
| `initial_lower` | 0.001 | 10.0 | `min(responses) * 0.1` |
| `initial_upper` | 3.784 | 90.0 | `max(responses) * 1.1` |
| `initial_ec50` | 108.0 | 0.4 | `median(concentrations)` |
| `initial_hill_coefficient` | 1.515 | -0.3 | N/A |
| `maxfev` | 3,000,000 | 3,000,000 | N/A |
| `zero_replacement` | 1e-24 | 1e-24 | N/A |
| `curve_direction` | Auto-detect (tries both) | Auto-detect (tries both) | N/A |

## Processing-level options

The parameters above control **curve fitting** (e.g. `hill_fitting.fit_hill_curve`, or `**fit_params` passed to `sp.process`). A separate **processing** option affects when pre-calc Hill params are overwritten by fitted values:

- **`allow_overwrite_hill_coefficients`**: When your CSV has **both** raw DATA/CONC **and** pre-calc (AC50, Upper, Lower, Hill_Slope, r2), sprime fits from raw and overwrites pre-calc. By default it **raises** unless you set `allow_overwrite_hill_coefficients=True`. When you allow overwrite, sprime **logs a warning** (console + report) each time pre-calc is overwritten.

See [Basic Usage Guide](basic_usage_guide.md) § "Processing option: allow_overwrite_hill_coefficients" for details.

## References

For more information on the underlying optimization algorithm and 4PL model, see:
- [Understanding 4PL Dose-Response Curves](../README_4PL_Dose_Response.md) - User-friendly explanation of 4PL model
- [Background and Concepts](../background_and_concepts.md) - qHTS background and terminology
- [Basic Usage Guide](basic_usage_guide.md) - How to use sprime with your data

