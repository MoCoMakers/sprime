# Terminology reference (sprime)

Plain-language definitions for words and ideas that usually **aren't taught in elementary school**, plus terms that are **specific to this library and the S' workflow**. If a word is everyday English, it is usually omitted here.

**Related docs:** [Background and Concepts](../background/background_and_concepts.md) * [S' derivation pipeline](../background/s_prime_derivation_pipeline.md) * [4PL / Hill curves](../background/README_4PL_Dose_Response.md) * [Basic usage guide](basic_usage_guide.md)

---

## Screening and biology (general)

| Term | Simple explanation |
|------|--------------------|
| **Preclinical** | Research done *before* human trials (often in dishes or animals) to see if an idea is safe or effective enough to test further. |
| **Drug discovery** | The process of finding and testing molecules that might become medicines. |
| **Quantitative high-throughput screening (qHTS)** | Testing **many** compounds at once, at **several doses**, and recording **numbers** (not just "yes/no") for how cells or targets respond. |
| **Assay** | A defined lab recipe: what you measure, how you measure it, under what conditions (time, plate, reader, etc.). |
| **Cell line** | Cells grown in the lab that came from one source and are treated as one "type" for experiments (e.g. a tumor line vs a non-tumor line). |
| **Compound** | Here: a small-molecule drug candidate (or similar) being tested--not necessarily an approved medicine yet. |
| **Dose / concentration** | **How much** of the compound the cells see (often in **molar** units like micromolar, **uM**). |
| **Dose-response** | How the measured **response** changes as **concentration** increases--often an S-shaped curve. |
| **Readout** | The number the instrument reports (fluorescence, luminescence, "% activity," viability, etc.). |
| **DMSO** | A common **solvent** used to dissolve compounds on screening plates. Wells with "only DMSO" (no drug) act as a **vehicle control**. |
| **Vehicle control** | The "no compound, same solvent" condition used to compare compound wells. In prose we often say **DMSO** even if another vehicle was used. |
| **Potency** | **How little** drug you need to get an effect (often summarized by **EC50**). |
| **Efficacy** | **How strong** the effect can get at high dose (related to curve **plateau / asymptotes**). |
| **In vitro** | In a dish or tube (outside a living body). |

---

## Math and curve fitting (general)

| Term | Simple explanation |
|------|--------------------|
| **Sigmoid / S-shaped curve** | A curve that is flat at low dose, rises (or falls) steeply in the middle, then levels off again. |
| **Four-parameter logistic (4PL)** | A standard **S-shaped** model with four parameters. Generic texts often say "bottom/top"; in sprime the Hill implementation uses **zero_asymptote** and **inf_asymptote** (dose-axis **left** vs **right** limits), **EC50**, and **steepness_coefficient**. |
| **Hill curve / Hill equation** | The specific 4PL-style formula sprime fits to dose-response points (see [4PL doc](../background/README_4PL_Dose_Response.md)). |
| **EC50 / AC50** | Concentration where the response is **halfway** between the two asymptotic limits (for a standard sigmoid). "AC50" and "EC50" are used interchangeably in many files. |
| **Asymptote** | The value the curve **approaches** as concentration -> 0 or at saturating dose (flat regions of the sigmoid). |
| **Zero asymptote** | Response in the limit **concentration -> 0** -- **left** side of the dose axis in sprime's Hill form. Not necessarily the numerically smaller response (curves may go up or down with dose). |
| **Inf asymptote** | Response in the limit **saturating concentration** -- **right** side of the dose axis. |
| **Hill slope / steepness** | How **sharp** the transition is between the two asymptotic limits (in code: often `steepness_coefficient`). |
| **R^2 (R-squared)** | A 0-1 score for how well the curve **matches** the points (1 = perfect fit; low values warn of noisy or non-sigmoidal data). |
| **Log scale / linear scale** | **Linear x** here means concentrations are used in a form suitable for the implemented fit (see fitting docs); "log" in file headers often refers to **log10(concentration)** for display. |

---

## S' and comparison metrics (sprime-specific)

| Term | Simple explanation |
|------|--------------------|
| **S' (S prime)** | One number that **summarizes** a dose-response curve using **potency and efficacy together** (from asymptotes and EC50). **Higher S'** usually means a **stronger** combined profile in this framework. Formula: **asinh((zero_asymptote - inf_asymptote) / EC50)** ("Generation 2"; see [Background](../background/background_and_concepts.md#the-s-s-prime-metric)). Use **fitted** or **file** asymptotes — do not assume values such as ``inf_asymptote = 0``. If a table gives only **span** (A-D) and **EC50**, use **asinh(span/EC50)**. |
| **Generation 2 (S')** | The current sprime definition of the score, evolved from an earlier **S** metric (see README citation). |
| **Delta S' (ΔS′)** | **Difference** between S' in a **reference** cell line and a **test** cell line: `S'(reference) - S'(test)`. Used to discuss **selectivity** between contexts. |
| **Profile** | One compound + cell line (+ assay context) curve--one row's worth of dose-response or pre-calculated parameters. |

---

## Two different uses of "normalized" (important)

In this project **"normalized" does not mean one thing**. Always check which step is meant:

| Phrase | Meaning |
|--------|--------|
| **Control-response normalization** (DMSO / vehicle) | Divide raw responses by a **control readout** for the same experimental context (**`Control_Response`**). Puts curves on a **relative** scale and reduces **plate-to-plate** drift. |
| **Asymptote-normalized / max scale** | After (optional) control step, **rescale** so the **largest** value is **1** (others ≤ 1), then often **x100** before fitting. **Different** from the control ratio. |

See [S' derivation pipeline Sec.3-4](../background/s_prime_derivation_pipeline.md).

---

## Pipeline names and code (sprime-specific)

| Term | Meaning |
|------|--------|
| **`sprime.response_pipeline`** | Python module that implements **shared math** for control ratios, optional **normalize_to_max_value** (then x100), and the **x100** convention--aligned with the **variation reference** spreadsheet columns (test/control ratio, normalize-to-1, x100, nonnormalized x100). |
| **`ratios_to_control`** | `response / Control_Response` (each point divided by the vehicle control readout). |
| **`normalize_to_max_value`** | After control ratios, scale so the **maximum** value is **1** and all other values are proportionally **≤ 1** (then typically ×100). Often the max is at low dose on descending curves, matching spreadsheet *normalise to 1*. |
| **`scale_responses`** | Multiply by **`S_PRIME_RESPONSE_SCALE_FACTOR`** (100) after max normalization--**convention** for readable magnitudes, not a biological constant by itself. |
| **`pipeline_response_scale`** | Control ratio -> **x100**, **without** max normalization (**response-scale** / "nonnormalized" path in the reference sheet). |
| **`pipeline_asymptote_normalized`** | Control ratio -> **normalize_to_max_value** -> **x100** (**asymptote-normalized** path). |
| **Raw path** | You have **concentrations + responses** -> library **fits** Hill parameters -> then S'. |
| **Pre-calculated path** | File already has **EC50 + asymptotes** (and maybe slope, R^2) -> S' computed **without** refitting. |
| **`skip_control_response_normalization`** | Import-time flag: if `True`, **skip** dividing by **`Control_Response`**. Default **`False`** (use control when provided). **Discouraged** for true raw plate readouts unless data are already relative. |
| **`allow_overwrite_precalc_params`** | If a row has **both** raw points **and** pre-calculated Hill params, sprime normally **refuses** to overwrite; set **`True`** to refit from raw and replace pre-calculated params (warnings are logged). |

---

## Data files and columns (sprime-specific)

| Term | Meaning |
|------|--------|
| **Path A** | **Raw** dose-response data (`Data*`/`Conc*` columns, or list format). |
| **Path B** | **Pre-calculated** Hill parameters (`AC50`/`ec50`, asymptotes, etc.). |
| **`values_as="columns"` / `"list"`** | Raw layout: one column per concentration/response vs comma-separated pairs in **`Responses`** / **`Concentrations`**. |
| **`Concentration_Units`** | Required for Path A; values are converted internally to **uM** (see supported synonyms in README / basic guide). |
| **`Control_Response`** | Column for the **vehicle control** readout used with raw curves when applying control-response normalization. |
| **Reserved vs pass-through columns** | Some headers are **reserved** for the pipeline (e.g. `Cell_Line`, `Compound_ID`). Others (e.g. MOA) are **metadata** carried through to exports when configured. |
| **Variation reference CSV** | `tests/fixtures/SPrime_variation_reference.csv` -- committed **test fixture** with **parallel** normalized vs response-scale columns; used to validate `response_pipeline`. **Human-maintained**; small rounding differences vs strict float pipelines can occur. |

---

## Optional identifiers and metadata (abbreviations)

| Term | Meaning |
|------|--------|
| **SMILES** | A text line describing **chemical structure** for tools/chemists. |
| **PubChem SID** | PubChem **substance** identifier. |
| **NCGC / NCATS** | Common sources of screening compound IDs (context-dependent). |
| **MOA** | Mechanism of action (how a drug is thought to work)--often a **label**, not used to compute S'. |

---

## Quick "where do I read more?"

| Question | Doc |
|----------|-----|
| Big picture on qHTS and S' | [Background and Concepts](../background/background_and_concepts.md) |
| Branch diagram: raw vs precalc, DMSO, max norm, x100 | [S' derivation pipeline](../background/s_prime_derivation_pipeline.md) |
| Hill equation details | [README_4PL_Dose_Response.md](../background/README_4PL_Dose_Response.md) |
| CSV formats, `load`/`process`, delta S' | [Basic usage guide](basic_usage_guide.md) |
| Fitter knobs | [Hill curve fitting configuration](hill_curve_fitting_configuration.md) |

---

*This glossary is for readers; the **code and tests** remain the implementation reference.*
