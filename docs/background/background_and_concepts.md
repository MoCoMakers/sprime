# Background and Concepts

This document provides background information on quantitative high-throughput screening (qHTS), assays, and the S' metric for users of the sprime library.

## Table of Contents

1. [Quantitative High-Throughput Screening (qHTS)](#quantitative-high-throughput-screening-qhts)
2. [Understanding Assays](#understanding-assays)
3. [Key Concepts and Terminology](#key-concepts-and-terminology)
4. [S' preprocessing branches (overview)](#s-preprocessing-branches-overview)
5. [The S' (S Prime) Metric](#the-s-s-prime-metric)
6. [Delta S' for Comparative Analysis](#delta-s-for-comparative-analysis)

## Quantitative High-Throughput Screening (qHTS)

Quantitative High-Throughput Screening (qHTS) is a method used in preclinical drug discovery to systematically test large libraries of compounds against biological targets, typically cell lines. The goal is to identify compounds with therapeutic potential by measuring dose-response relationships across multiple concentrations.

### Key Characteristics

- **Systematic testing**: Large numbers of compounds tested in parallel
- **Multiple concentrations**: Each compound tested at various doses/concentrations
- **Quantitative measurements**: Responses measured numerically (e.g., viability, activity, luminescence)
- **Dose-response relationships**: Analysis focuses on how response changes with concentration

## Understanding Assays

### What is an Assay?

An **Assay** is a standardized biological procedure that defines:

- The experimental methodology and protocol
- The equipment and machinery used
- The measurement technique (e.g., luminescence, fluorescence, viability)
- The test conditions (temperature, time points, media, etc.)
- Quality control standards and acceptance criteria

### Assay Identification

Assays are typically identified in data tables through:

- **Screen ID** (e.g., `HTS002`) - a unique identifier for the screening campaign
- **Assay Name/ID column** - an optional identifier column in the data table
- **Readout type** (e.g., "activity", "viability", "luminescence") - the type of measurement

Multiple dose-response profiles can belong to the same assay if they were measured under identical experimental conditions.

### Time Profiles in Assays

Many assays include time-dependent measurements, where the response is measured at specific time points after compound exposure:

- **24Hr** - 24-hour time point
- **48Hr** - 48-hour time point  
- **4Day** - 4-day (96-hour) time point

Data at each time point represents a separate dose-response profile. The same compound tested against the same cell line at different time points may show different responses, reflecting how the compound's effects evolve over time.

**Example:** A compound might show minimal effect at 24Hr but strong cytotoxicity at 4Day, indicating delayed response kinetics.

### Assay Standardization Importance

Within a single assay, all dose-response profiles are measured under identical conditions:
- Same time point (if time-course study)
- Same equipment and readout method
- Same control conditions
- Same quality control standards

This standardization allows:
- Direct comparison between compounds tested in the same assay
- Reliable delta S' calculations (comparing reference vs test cell lines)
- Reproducible analysis and interpretation

When combining data from multiple assays, researchers must be cautious as differences in methodology can introduce variability that affects comparisons.

## Key Concepts and Terminology

For a **compact glossary** (including sprime-specific names like *asymptote-normalized* vs *response-scale*), see the **[Terminology reference](usage/terminology_reference.md)** in `docs/usage/`.

### Cell Line

A **cell line** is a clone or sampled tissue - frequently a TERT-mutated sample from a tissue biopsy - thought to be genetically identical within the assay. Cell lines represent biological systems being tested (e.g., tumor cells, normal cells, specific tissue types).

**Examples:**
- `ipNF96.11C` - A neurofibroma cell line
- `LS513_LARGE_INTESTINE` - A large intestine cancer cell line (CCLE format)
- `ipnf05.5 mc` - A reference cell line (e.g., non-tumor tissue)

### Compound

A **compound** is generally a small molecule drug being investigated for therapeutic effects. In qHTS context, compounds are tested against cell lines to measure their biological activity.

**Common identifiers:**
- Compound Name (e.g., "Trifluoperazine")
- Compound_ID (e.g., "NCGC00013226-15")
- NCGC ID - National Center for Advancing Translational Sciences (NCATS) identifier
- Broad ID - Broad Institute identifier
- PubChem SID - PubChem substance identifier
- SMILES - Chemical structure notation

### Hill Curve (4-Parameter Logistic Model)

A **Hill curve** is a four parameter logistic model for pharmaceutical compounds. The X-axis represents concentration, and the Y-axis represents some form of response. The response is often % response relative to DMSO control, but it may be direct luminescence response in a bioluminescence assay (of cell survival/replication).

The Hill equation models the sigmoidal (S-shaped) relationship between concentration and response:

```
y = D + (A - D) / (1 + (x/C)^n)
```

Where:
- **A** = **zero_asymptote** (response as concentration -> 0; left side of dose axis)
- **D** = **inf_asymptote** (response at saturating concentration; right side of dose axis)
- **C** = EC50 (half-maximal concentration)
- **n** = Hill coefficient (slope/steepness)

### R-squared (R^2)

**R-squared (R^2)** is a goodness of fit statistic historically used to assess how well dose-response data fits the Hill equation model. R^2 ranges from 0 to 1, where 1 = perfect fit. Higher R^2 indicates better curve quality.

Historically used to filter and validate dose-response profiles before analysis:

- **R^2 = 1.0** - Perfect fit
- **R^2 > 0.9** - Excellent fit (typically considered high quality)
- **R^2 0.7-0.9** - Good fit
- **R^2 < 0.7** - Poor fit (may indicate experimental issues or non-sigmoidal response)

Poor fit (low R^2) can indicate:
- Insufficient data points
- Data not following sigmoidal dose-response pattern
- Experimental artifacts or errors
- Saturation issues (concentrations not spanning full response range)

### Data Structure

#### Core Dose-Response Data

Each dose-response profile contains:

1. **Raw Response Data** - Measured responses at various concentrations
   - Column naming: `DATA0`, `DATA1`, `DATA2`, ..., `DATAN`
   - Values can be: % response relative to control, absolute luminescence/fluorescence, viability percentage, etc.
   
2. **Concentration Data** - Compound concentrations tested
   - Column naming: `CONC0`, `CONC1`, `CONC2`, ..., `CONCN`
   - Units: Typically microMolar (uM), but may vary
   - Format: Scientific notation common (e.g., `1.30E-09` = 1.30 x 10-9 M)

3. **Fitted Hill Curve Parameters** (may be pre-calculated or calculated from raw data)
   - **AC50** or **EC50** - Half-maximal activity/effect concentration
   - **Inf_asymptote** (or legacy Upper/Infinity) - `inf_asymptote` (saturating dose, right of dose axis)
   - **Zero_asymptote** (or legacy Lower/Zero) - `zero_asymptote` (concentration -> 0, left of dose axis)
   - **Hill** (slope) - Hill coefficient (steepness of curve)
   - **R^2** (r2) - Coefficient of determination (goodness of fit)

#### Additional Metadata

Common fields across different data sources:

- **Compound Information**: Name, Compound_ID, NCGCID (optional pass-through), SMILES, PubChem SID
- **Cell Line Information**: Name (e.g., `ipNF96.11C`), Ref ID (e.g., `ACH-000007`, `depmap_id`)
- **Biological Context**: Target(s), MoA (Mechanism of Action), Disease Area, Indication, Phase
- **Quality Metrics**: R^2, AUC (Area under the curve), MAXR (Maximum response), Classification fields (CCLASS, CCLASS2)

## S' preprocessing branches (overview)

The **same S' formula** (see [below](#the-s-s-prime-metric)) can be fed by different **upstream pipelines**. Understanding these branches avoids confusion when people say "normalized" in screening workflows.

| Branch | Short description |
|--------|-------------------|
| **Pre-calculated vs raw** | Some files already contain EC50 and asymptotes; others provide only concentration-response series. sprime **prefers raw** data when possible so curves can be **fit** and QC'd in one place; pre-calculated parameters remain supported when raw data are unavailable. |
| **DMSO / vehicle step (first "normalized")** | For **raw** readouts, responses are often expressed as **ratio to a reference vehicle control** per sample (\(r/d\)); in CSV this is **`Control_Response`** (*DMSO* is biology shorthand). **Best practice** is a control matched to the **same experimental context** (e.g. per plate). Use **`skip_control_response_normalization=True`** at **import** when responses are already on the analysis scale (empty **`Control_Response`** allowed); **\(r/d\)** math is **`response_pipeline`**, not `load`--see [S' derivation pipeline Sec.3.4](s_prime_derivation_pipeline.md). |
| **Baseline / asymptote step (second "normalized")** | After the DMSO step (or skip), you may **scale ratios so the maximum is 1** (then ×100; **asymptote-normalized** path--recommended for cross-dataset magnitude correction when the low-dose peak is the max) or **omit** that step (**response-scale** path--e.g. historical alignment with published NF1 / PNF1 analyses). |
| **x100 convention** | Both major raw variants typically apply a **fixed x100** multiplier for **numerical convenience**; it is **not** biologically special by itself. |

**Full detail (glossary, order of operations, implementation notes):** [S' derivation pipeline](s_prime_derivation_pipeline.md).

## The S' (S Prime) Metric

S' is a single value score that summarizes a drug's dose-response curve relative to a single cell line and assay.

### Calculation

S' is calculated from Hill curve parameters using the formula:

**S' = asinh((Zero_asymptote - Inf_asymptote) / EC50)**

This matches the **named** Hill parameters (same order as CSV **Zero**/**Lower** vs **Inf**/**Upper**). It is **not** the same as :math:`(\mathrm{Inf}-\mathrm{Zero})/\mathrm{EC50}` inside ``asinh``--the 4PL fit still stores ``inf_asymptote`` and ``zero_asymptote`` in Hill form; S' uses their **difference in this order** so the sign tracks whether response **drops** (typical inhibition: Zero > Inf => positive S') or **rises** along dose.

This is equivalent to:

$$S' = \ln\left(\frac{\text{Zero\_asymptote} - \text{Inf\_asymptote}}{\text{EC50}} + \sqrt{\left(\frac{\text{Zero\_asymptote} - \text{Inf\_asymptote}}{\text{EC50}}\right)^2 + 1}\right)$$

Where:
- **Zero asymptote** = `zero_asymptote` (concentration -> 0, left of dose axis)
- **Inf asymptote** = `inf_asymptote` (saturating dose, right of dose axis)
- **EC50** = Half-maximal effect concentration (C parameter)

### What S' Represents

- **Larger magnitude |S'|** reflects a stronger combined potency/efficacy signal (given consistent curve orientation and preprocessing)
- **Sign** depends on asymptote ordering (inhibition vs activation shapes); compare compounds under the same assay convention
- **S' combines potency and efficacy** into a single metric, making it useful for ranking compounds

### Dose-Response Curve Requirements

For S' to be meaningful, dose-response curves must have:
- Sufficient data points across a relevant range of concentrations
- Data modeled with a four-parameter logistic curve (Hill Equation)
- Resulting in reliable estimates of Zero asymptote, Inf asymptote, and EC50

### Generation 2 Methodology

This library implements **Generation 2** of the S' methodology, which evolved from the original S metric described in:

> Zamora PO, Altay G, Santamaria U, Dwarshuis N, Donthi H, Moon CI, Bakalar D, Zamora M. Drug Responses in Plexiform Neurofibroma Type I (PNF1) Cell Lines Using High-Throughput Data and Combined Effectiveness and Potency. *Cancers (Basel)*. 2023 Dec 12;15(24):5811. doi: [10.3390/cancers15245811](https://doi.org/10.3390/cancers15245811). PMID: 38136356; PMCID: PMC10742026.

New citations introducing Generation 2 (S') are forthcoming.

## Delta S' for Comparative Analysis

**Delta S'** (ΔS') enables quantitative comparison of drug responses between different cell lines within a single assay. It is calculated as:

**ΔS' = S'(reference cell line) - S'(test cell line)**

### Use Cases

This metric allows researchers to:

- **Compare drug effects across cell lines**: Identify compounds that show differential responses between reference (e.g., non-tumor tissue) and test cell lines (e.g., tumor cell lines)
- **Rank compounds by selectivity**: Compounds with higher (more negative) delta S' values indicate greater differential effects in the test cell line relative to the reference
- **Prioritize drug candidates**: Within a single assay, compounds can be ranked by their delta S' values, helping identify the most promising candidates for further investigation

### Example Calculation

If a compound has:
- S' = 2.5 for a reference cell line (normal tissue)
- S' = 4.0 for a test cell line (tumor tissue)
- Then ΔS' = 2.5 - 4.0 = -1.5

This indicates the compound is more effective (higher S') in the tumor cell line.

### Interpreting Delta S'

- **Negative delta S' values**: More effective in test cell line (higher S' in test than reference)
- **Positive delta S' values**: More effective in reference cell line (higher S' in reference than test)
- **More negative = more selective**: Compounds with more negative delta S' values are typically more selective for the test (e.g., tumor) cell line and may be prioritized for further investigation

A ranking system based on delta S' values allows researchers to systematically identify compounds with the greatest therapeutic potential within their screening assay.

## Further Reading

- [S' derivation pipeline](s_prime_derivation_pipeline.md) - **Canonical** narrative: pre-calculated vs raw, DMSO-relative ratios, two meanings of "normalized," asymptote-normalized vs response-scale, x100 convention
- [Basic Usage Guide](../usage/basic_usage_guide.md) - How to use sprime with your data
- [Understanding 4PL Dose-Response Curves](README_4PL_Dose_Response.md) - Detailed explanation of the Hill equation model
- [Hill Curve Fitting Configuration](../usage/hill_curve_fitting_configuration.md) - Technical details on curve fitting parameters





