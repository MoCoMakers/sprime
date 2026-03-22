# Understanding 4PL Dose-Response Curves

## Introduction

This guide explains the **Four Parameter Logistic (4PL) model** and dose-response curves for researchers and biologists using the sprime library. You don't need to be a programmer to understand these concepts--this guide focuses on the biological and statistical meaning of dose-response analysis.

## What is a Dose-Response Curve?

A **dose-response curve** (also called a concentration-response curve) shows how a biological system responds to different concentrations of a compound (drug, chemical, etc.). The curve typically has a sigmoidal (S-shaped) appearance:

- **Low concentrations**: Minimal or baseline response
- **Intermediate concentrations**: Steep increase (or decrease) in response
- **High concentrations**: Response plateaus at maximum (or minimum) level

### Why Are Dose-Response Curves Important?

Dose-response curves are fundamental to:
- **Drug discovery**: Understanding how effective a compound is at different doses
- **Toxicology**: Determining safe concentration ranges
- **Pharmacology**: Characterizing drug potency and efficacy
- **High-throughput screening**: Ranking compounds by their biological activity

## The Four Parameter Logistic (4PL) Model

The 4PL model is a mathematical equation that describes sigmoidal dose-response relationships. It's also known as the **Hill equation**, named after Archibald Hill who developed it to describe oxygen binding to hemoglobin.

### The Mathematical Formula

**sprime** implements the linear-x Hill / 4PL in the same form as `hill_equation` in code (parameter names match `HillCurveParams`):

```
y = inf_asymptote + (zero_asymptote - inf_asymptote) / (1 + (x / ec50) ** steepness_coefficient)
```

Equivalent display (ASCII-only identifiers in the fraction):

$$y = \mathrm{inf\_asymptote} + \frac{\mathrm{zero\_asymptote} - \mathrm{inf\_asymptote}}{1 + \left(\frac{x}{\mathrm{EC50}}\right)^{n}}$$

Where:
- **y** = Response (what you measure: % viability, fluorescence, activity, etc.)
- **x** = Concentration (same units as EC50)
- **zero_asymptote** = Response approached as **concentration approaches zero**
- **inf_asymptote** = Response approached as **concentration becomes large** (saturating dose region)
- **n** = `steepness_coefficient` in code (Hill exponent in `(x/EC50)^n`; not the same number as log-x tools often call "Hill slope")
- **EC50** = Half-maximal concentration (`ec50` in code)

Either asymptote may be numerically larger depending on whether the curve rises or falls with dose.

Many textbooks write the same shape as **y = D + (A - D) / (...)** with **A** = limit at low **x** and **D** at high **x**; that matches **zero_asymptote** = A and **inf_asymptote** = D for the usual sign of **n** in each curve direction.

### The Four Parameters Explained

#### 1. Zero asymptote (`zero_asymptote`) -- limit as concentration approaches zero

**What it means:**
- The response value approached when concentration is very low (vehicle / lowest doses in the series)
- Not named for being numerically "low" or "high": for a **decreasing** dose-response curve, `zero_asymptote` is often **larger** than `inf_asymptote`; for an **increasing** curve, the opposite may hold

**Example:**
- At the vehicle / lowest-dose end of the series, the fitted plateau might be ~100 on a % scale while the high-dose plateau is ~5 -- here **100** is still `zero_asymptote`

#### 2. Inf asymptote (`inf_asymptote`) -- limit at saturating concentration

**What it means:**
- The response value approached when concentration is very high (saturating dose region)
- Not "upper" or "lower" in value -- only **which dose extreme** the parameter describes

**Example:**
- Residual signal at the highest tested concentrations (e.g. ~0.85 on a x100 vs control scale) is `inf_asymptote` even when that number is small

#### 3. EC50 (`ec50` in code) - Half-Maximal Concentration

**What it means:**
- The concentration at which the response is halfway between the zero asymptote and the inf asymptote
- **EC50** = Effective Concentration 50% (for activation)
- **IC50** = Inhibitory Concentration 50% (for inhibition)
- A key measure of **potency**: lower EC50 = more potent compound

**Example:**
- If EC50 = 10 uM, this means at 10 micromolar concentration, the response is at 50% of the maximum range
- A compound with EC50 = 1 uM is more potent than one with EC50 = 100 uM

**Why it matters:**
- EC50 is the most commonly reported value in drug screening
- Allows comparison of compound potency across different assays
- Lower EC50 values indicate compounds that work at lower concentrations

#### 4. Steepness coefficient (`steepness_coefficient`, exponent **n**) - Steepness

**What it means:**
- Describes how steeply the curve transitions between `zero_asymptote` and `inf_asymptote`
- In other tools this is often loosely called "slope" or "Hill slope" (not the same number as log-x Hill slope)
- Higher absolute values = steeper curve = more cooperative response
- In **sprime** (linear-x 4PL, concentration on a linear scale in `(x/EC50)^n`): positive values = increasing curve (activation); negative values = decreasing curve (inhibition). Tools that use a **log-x 4PL** (log10 concentration in the exponent) may report a positive Hill slope for the same biological curve; the models are not interchangeable without conversion. See **[Linear-x vs log-x 4PL (Hill slope)](#linear-x-vs-log-x-4pl-hill-slope)** below.

**Example:**
- n = 1.5: moderate steepness, typical for many biological systems
- n = 3.0: very steep curve, suggests cooperative binding
- n = 0.5: shallow curve, suggests less cooperative response

**Biological interpretation:**
- Exponent n about 1: Simple binding (one molecule per target)
- n > 1: Cooperative binding (multiple molecules work together)
- n < 1: Negative cooperativity or multiple binding sites

## Linear-x vs log-x 4PL (Hill slope)

**sprime** fits a **linear-x** four-parameter logistic: concentration *x* (linear scale, same units as EC50) enters the Hill form as `(x/EC50)^n`. The exponent *n* is **`steepness_coefficient`** in code (the same mathematical role as the classical **Hill coefficient** *n* in many linear-x texts).

Many graphing and pharmacology tools use a **log-x** 4PL instead (e.g. base-10 log of concentration in the exponent). They often report a "**Hill slope**" that is **not numerically the same** as sprime's *n* for the same biological curve. Do not mix slope values across parameterizations without conversion. For implementation notes, see the **`sprime.hill_fitting`** module docstring.

## Visualizing Dose-Response Curves

### Typical Curve Shapes

**Increasing Curve (Activation):**
```
Response
  100% |                    /------- inf_asymptote
       |                  /
   50% |                /  <- EC50 point
       |              /
    0% |-------------/
       |         zero_asymptote
       +------------+------------+------------
         0.01       EC50         100
              Concentration (uM)
```

**Decreasing Curve (Inhibition):**
```
Response
  100% |-------------\  zero_asymptote
       |              \
   50% |                \  <- IC50 point
       |                  \
    0% |                    \------- inf_asymptote
       +------------+------------+------------
         0.01       IC50         100
              Concentration (uM)
```

## R-Squared: How Well Does the Curve Fit?

**R^2 (R-squared)** is a statistical measure that tells you how well your experimental data points match the fitted 4PL curve.

- **R^2 = 1.0**: Perfect fit (all data points lie exactly on the curve)
- **R^2 = 0.9**: Very good fit (90% of variance explained)
- **R^2 = 0.7**: Moderate fit (70% of variance explained)
- **R^2 < 0.5**: Poor fit (data may not follow sigmoidal pattern)

**What affects R^2?**
- **Data quality**: Noisy or scattered data, so lower R^2
- **Data range**: Need concentrations spanning the full response range
- **Outliers**: Single bad data points can reduce R^2
- **Biological variability**: Some assays are inherently more variable

**When is R^2 acceptable?**
- **R^2 > 0.9**: Excellent, publishable quality
- **R^2 > 0.8**: Good, reliable for most purposes
- **R^2 > 0.7**: Acceptable, but may want to investigate data quality
- **R^2 < 0.7**: Consider re-running assay or checking for issues

## How sprime Uses 4PL Fitting

When you use sprime to analyze your screening data:

1. **Input**: You provide raw dose-response data (concentrations and responses)

2. **Fitting**: sprime automatically fits a 4PL curve to your data using mathematical optimization

3. **Output**: You get:
   - **EC50**: Potency value
   - **zero_asymptote** and **inf_asymptote**: response limits at low and high concentration
   - **steepness_coefficient** (exponent n): curve steepness
   - **R^2**: fit quality

4. **S' Calculation**: sprime then calculates the **S' (S prime)** metric from these parameters, which provides a single value summarizing the entire dose-response profile

## Common Questions

### Why use 4PL instead of simpler models?

The 4PL model captures the full sigmoidal shape of biological dose-response relationships. Simpler models (like linear) don't account for:
- Baseline responses at low concentrations
- Maximum achievable responses at high concentrations
- The S-shaped transition between these extremes

### What if my data doesn't look sigmoidal?

If your data doesn't follow a sigmoidal pattern, the 4PL fit may have low R^2. Possible reasons:
- **Insufficient concentration range**: Need data from baseline to plateau
- **Non-sigmoidal biology**: Some responses are linear or have different shapes
- **Data quality issues**: Outliers, experimental errors, or assay problems

### Can I compare EC50 values across different assays?

Yes, but be cautious:
- EC50 values are comparable within the same assay type
- Different readouts (viability vs. activity) may have different scales
- Always consider the biological context and assay conditions

### What's a "good" EC50 value?

It depends on your context:
- **Drug discovery**: Often looking for EC50 < 1 uM (nanomolar range is ideal)
- **Toxicology**: May be interested in higher concentrations
- **Research**: Depends on your specific question

The important thing is comparing EC50 values **within your screening assay** to rank compounds.

## Real-World Example

Imagine you're screening compounds for anti-cancer activity:

**Compound A:**
- EC50 = 0.5 uM (very potent!)
- R^2 = 0.95 (excellent fit)
- inf_asymptote = 90% (strong effect at high dose)

**Compound B:**
- EC50 = 50 uM (less potent)
- R^2 = 0.88 (good fit)
- inf_asymptote = 75% (moderate effect at high dose)

**Interpretation:**
- Compound A is **100x more potent** than Compound B (0.5 vs 50 uM)
- Compound A also has a **stronger maximum effect** (90% vs 75%)
- Both have good curve fits (R^2 > 0.85)
- **Conclusion**: Compound A is the better candidate

## Further Reading

For more information:
- **Configuration options**: See [Hill Curve Fitting Configuration](usage/hill_curve_fitting_configuration.md)
- **Background and concepts**: See [Background and Concepts](background_and_concepts.md) for qHTS and assay terminology
- **Theoretical foundations**: See references below

## References

1. **Hill, A. V.** (1910). The possible effects of the aggregation of the molecules of haemoglobin on its dissociation curves. *J. Physiol.*, 40, iv-vii.

2. **De Lean, A., Munson, P. J., & Rodbard, D.** (1978). Simultaneous analysis of families of sigmoidal curves: application to bioassay, radioligand assay, and physiological dose-response curves. *Am. J. Physiol.*, 235(2), E97-E102.

3. **Gottschalk, P. G., & Dunn, J. R.** (2005). The five-parameter logistic: a characterization and comparison with the four-parameter logistic. *Anal. Biochem.*, 343(1), 54-65.

4. **Cardillo, G.** (2012). Four parameters logistic regression - There and back again. MATLAB Central File Exchange. Retrieved from https://it.mathworks.com/matlabcentral/fileexchange/38122

5. **AAT Bioquest, Inc.** (2024). Quest Graph(TM) Four Parameter Logistic (4PL) Curve Calculator. Retrieved from https://www.aatbio.com/tools/four-parameter-logistic-4pl-curve-regression-online-calculator

