# S' derivation pipeline: branches, controls, and normalization (conceptual reference)

This document is the **canonical conceptual description** of how sprime's S' workflow is intended to branch: what inputs are possible, what **best practices** we recommend, where the word **"normalized"** means different things, and how the **x100 scale** convention fits in. It complements [Background and Concepts](background_and_concepts.md) and the [Basic Usage Guide](../usage/basic_usage_guide.md).

**Related materials**

- Variation reference (committed test fixture, parallel calculation paths): [tests/fixtures/SPrime_variation_reference.csv](../../tests/fixtures/SPrime_variation_reference.csv).

---

## 1. Why this document exists

Screening data arrive in **different shapes**. The same mathematical S' formula can be applied at the end, but **how** you obtain Hill parameters (zero asymptote, inf asymptote, EC50) depends on:

1. Whether the file already contains **pre-calculated** curve parameters or only **raw** dose-response points.
2. Whether raw responses are adjusted with a **reference DMSO (or vehicle) control** per sample.
3. Whether, after that step, you apply a **second normalization** that pegs the **baseline (low-concentration) response** to a reference level for **magnitude-corrected** comparison across datasets--or you keep the **historical response-scale** path.

sprime **encourages** using **raw** concentration-response data when possible so that the library (or your preprocessing) can **estimate** EC50 and asymptotes by fitting. Many legacy files, however, only ship **pre-fitted** parameters; those must remain supported.

---

## 2. Branch A: Pre-calculated parameters vs raw dose-response data

### 2.1 Pre-calculated path

Some datasets provide, per compound-cell line (and time point, if applicable), numeric fields such as:

- **EC50** / AC50 (half-maximal concentration)
- **Zero_asymptote** (response as concentration -> 0; left side of the dose axis -- not necessarily the numerically smaller response)
- **Inf_asymptote** (response at saturating concentration; right side of the dose axis)
- Optionally **Hill slope** / steepness, **r^2**, etc.

When these are **trusted** and no raw `DATA*` / `CONC*` (or list columns) are available, S' can be computed **directly** from the formula (see [Sec.7](#7-the-s-formula-scale-and-generation-2)) without fitting. This path does **not** by itself apply DMSO ratio or max normalization inside the library unless those effects are **already reflected** in the uploaded numbers.

**Limitation:** You cannot refit or QC the curve from raw points; you inherit whatever modeling choices the upstream pipeline used.

### 2.2 Raw path (preferred when possible)

When **raw** responses and concentrations are present, the recommended workflow is:

1. Apply **optional** DMSO-relative scaling ([Sec.3](#3-branch-b-dmso-reference-and-first-use-of-normalized)) unless explicitly skipped.
2. Apply **optional** max normalization ([Sec.4](#4-branch-c-second-normalization-baseline-pegging-vs-response-scale)) according to the chosen S' variant.
3. Fit a **linear-x 4PL / Hill** model to obtain zero asymptote, inf asymptote, EC50 (and steepness).
4. Apply the **x100** convention where required ([Sec.5](#5-the-100-scale-factor-convention)).
5. Compute S' from the fitted parameters ([Sec.7](#7-the-s-formula-scale-and-generation-2)).

**Encouragement:** Prefer raw data so that fitting, QC (e.g. R^2), and consistent modeling live in one place. **Not all studies** provide raw exports; the pre-calculated path remains first-class.

---

## 3. Branch B: DMSO reference and (first) use of "normalized"

### 3.1 Biological motivation

In plate-based screening, **absolute** readouts (luminescence, fluorescence counts, viability, etc.) drift between plates, instruments, and days. A standard mitigation is to express compound responses **relative to a vehicle control** on the **same experimental unit**--commonly **DMSO** (or the same solvent vehicle used for compound dilutions).

We use **"DMSO"** as shorthand for **that vehicle reference**; your assay might use a different control label, but the **role** is the same: a **reference raw response** to divide by (or subtract from, in other formulations--not sprime's default ratio form here).

### 3.2 Definition (DMSO-relative ratio)

For each sample (compound-cell line-plate-timepoint as your design requires), let:

- \(r\) = raw measured response at a given concentration (or the vector of responses before curve modeling).
- \(d\) = **reference DMSO / vehicle response** for that same experimental context.

A **DMSO-normalized ratio** (first sense of **"normalized"**) is conceptually:

\[
\text{ratio} = \frac{r}{d}
\]

(Exact implementation may use the same ratio for all concentrations in a curve after choosing one DMSO per curve; see study design.)

This step **controls for experimental variation** in overall signal level **across plates or batches** when DMSO is matched to the same context as the compound wells.

### 3.3 Best practice: granularity of DMSO

**Best practice** for screening assays is to attach a **reference DMSO value per sample** at a granularity that matches how variation occurs:

- **Ideal:** **Per plate** (or per replicate block)--DMSO measured alongside the same plate's compound wells.
- **Acceptable:** Per assay batch or study-wide default **only** when plate effects are negligible or the study design forces it.

The library's data model should allow **one control reference per profile** (or per row). The canonical CSV column name for that readout is **`Control_Response`** (*DMSO* in prose = vehicle control shorthand).

### 3.4 Import-time decision: strict raw load + `skip_control_response_normalization`

**What users should assume today (implementation, not only "target"):**

1. **Validation at load** -- For raw Path A, when **`skip_control_response_normalization=False`** (default), the CSV must include a **`Control_Response`** column and each raw row must carry a **non-empty, non-zero** numeric value. The loader **checks** this; it does **not** silently substitute a default.
2. **No automatic \(r/d\) inside `load` / `process` yet** -- sprime does **not** currently multiply or divide your **`DATA*`** / **`Responses`** by **`Control_Response`** when building profiles or fitting. Values passed to the Hill fit are **what you put in the file** (after unit conversion for concentrations only). To match reference spreadsheet columns such as *test/control ratio*, *x 100 scale*, or *Nonnormalized, x100 scale*, you must preprocess with **`sprime.response_pipeline`** (e.g. **`pipeline_response_scale`**, **`pipeline_asymptote_normalized`**) **before** fitting, or bake that scaling into the numbers upstream.
3. **When to set `skip_control_response_normalization=True`** -- Use when responses are **already** on your analysis scale (DMSO or other control already folded in). Then **`Control_Response`** may be **present but empty**; the loader skips the strict control check.

**Product direction (may evolve):** Optionally apply the \(r / \text{Control\_Response}\) step automatically at import and thread that intent through **`process()`** / exports so users do not have to call **`response_pipeline`** manually--see [Sec.8](#8-implementation-status-for-developers) and [Sec.10](#10-open-integration-questions-for-product--api-design).

Docstrings should use **DMSO** liberally as biology shorthand for the **control solution**, while the **column/API** name stays **`Control_Response`** / **`skip_control_response_normalization`**.

Skipping control-response normalization is **discouraged** for true raw plate readouts unless values are already on an appropriate relative scale or you accept **uncorrected plate-to-plate drift** in the strict-validation sense only (the numeric curve is still **not** auto-scaled by the library today).

> **No environment variables:** As a Python library, configuration is **method parameters and documented variants**, not env-based toggles.

---

## 4. Branch C: Second normalization -- max normalization vs response scale

After the DMSO step ([Sec.3](#3-branch-b-dmso-reference-and-first-use-of-normalized)) **or** after a deliberate skip ([Sec.3.4](#34-requirement-vs-opt-out-skip_dmso_normalization)), a **second** design choice affects magnitude and cross-dataset comparability. Here **"normalized"** means something **different** from Sec.3: it is **not** the DMSO ratio; it is **max normalization** (largest ratio → 1, then ×100) on the response axis.

### 4.1 Terminology: two meanings of "normalized"

| Term in docs | Meaning |
|--------------|--------|
| **DMSO-relative / vehicle normalization** | \(r / d\) -- controls **plate-to-plate** (or batch) **signal level** drift. |
| **Baseline / asymptote normalization** (second sense) | Rescale responses so the **low-concentration (baseline) behavior** is pegged to a reference (conceptually **1** on the ratio scale before optional x100), so **curve magnitude** is comparable across datasets. |

Always disambiguate which normalization you mean in prose (e.g. "DMSO-normalized ratio" vs "baseline-pegged / asymptote-normalized curve").

### 4.2 Recommended default: **asymptote-normalized** (baseline pegging)

**Intent:** At **effective zero compound concentration** (vehicle-only region), the response should align with a **consistent baseline** after DMSO ratio. In practice, noise can make the **observed** low-dose response deviate (e.g. **110%** of DMSO when it "should" be ~100%). **Baseline pegging** rescales responses so the **zero-asymptote region** is anchored (conceptually to **1** before downstream scaling), which:

- Mitigates **assay artifacts** and **baseline drift** in the **vertical** direction (magnitude).
- Improves **cross-dataset** comparability when the same S' variant is used.

**Implementation concept (after DMSO ratio):**

- **Normalize to max:** divide every ratio by **max(ratios)** so the largest value is **1** and all others are **≤ 1** (then typically multiply by **100**; see [Sec.5](#5-the-100-scale-factor-convention)). On typical **descending** qHTS curves the max is at the **low-dose / vehicle** end, aligning with spreadsheet *normalise to 1*; if a later point is numerically larger, the library scales to that peak.

We may refer to this in code or API as:

- **`asymptote_normalized`** -- aligns with tests and naming in this document.

**Synonyms for prose:** "max normalization," "normalize to 1" (sheet language when max is at low dose). Public API: **`normalize_to_max_value`** then ×100; pipeline flag **`asymptote_normalized`**.

### 4.3 **Response-scale** path (no max normalization)

**Intent:** Keep responses on the scale produced after the **control-response ratio** (and x100), **without** dividing by max first. Useful for matching **reference spreadsheet** non-normalized columns and for teaching **response-scale** workflows.

**Behavior:** After \(r/d\) when applicable, responses enter the fit **without** **normalize_to_max_value**. EC50 and asymptotes are estimated on that **response scale**.

**API name:** e.g. **`response_scale`** (matches tests such as `test_*_response_scale`).

**When to use:** Sensitivity analysis, parity with specific published tables, or when scientists explicitly want **no** max-normalization step. See **[demonstration.ipynb](../usage/demonstration.ipynb)** (route 3 / `response_scale`).

### 4.4 Order of operations (conceptual)

1. **Optional:** Control-response ratio \(r/d\) unless `skip_control_response_normalization=True`.
2. **Choose S' variant:**
   - **`asymptote_normalized`** -- **normalize_to_max_value** -> ×100 -> fit -> S'.
   - **`response_scale`** -- ×100 only (no max step) -> fit -> S'.
3. **Apply x100** where required by convention ([Sec.5](#5-the-100-scale-factor-convention)).
4. Compute S' from fitted (or pre-calculated) parameters ([Sec.7](#7-the-s-formula-scale-and-generation-2)).

---

## 5. The x100 scale factor (convention)

In **both** major raw pipelines (**asymptote-normalized** and **response-scale**), the reference workflow applies an **arbitrary multiplicative factor of 100** to responses (or to the ratio chain) before or during curve fitting **as documented for that pipeline**.

**Purpose:** Numerical convenience -- S' values become easier to report and compare (often modest magnitudes). **The factor has no special biological meaning by itself**; it is **required by convention** for that pipeline so that results stay comparable to reference tables and historical outputs.

Do **not** interpret "x100" as a percent label unless your upstream semantics already define it that way.

---

## 6. Relationship to the variation reference CSV

The fixture [SPrime_variation_reference.csv](../../tests/fixtures/SPrime_variation_reference.csv) illustrates parallel paths:

- **DMSO row** -- vehicle reference for ratio.
- **Normalized path** -- DMSO ratio -> **normalize_to_max_value** (max → 1) -> x100 -> 4PL -> S'.
- **Non-normalized path** -- DMSO ratio -> x100 **without** max normalization -> 4PL -> S'.

Align naming in code and docs with **Sec.3-Sec.5** above to avoid confusing the two "normalized" steps.

---

## 7. The S' formula, scale, and Generation 2

After zero asymptote, inf asymptote, and EC50 are known (fitted or pre-calculated), **Generation 2** S' uses:

\[
S' = \operatorname{asinh}\left(\frac{\text{Inf\_asymptote} - \text{Zero\_asymptote}}{\text{EC50}}\right)
\]

(or equivalent natural log form). **How** those asymptotes were produced--pre-calculated file, DMSO ratio only, or DMSO ratio + asymptote normalization--must be **tracked** for interpretation.

See also [Background and Concepts -- The S' metric](background_and_concepts.md#the-s-s-prime-metric) and the [README](../../README.md).

---

## 8. Implementation status (for developers)

- **`sprime.response_pipeline`** implements **test/control ratio**, **normalize_to_max_value**, **x100**, and the composed pipelines **`pipeline_asymptote_normalized`** / **`pipeline_response_scale`**, aligned with **`tests/fixtures/SPrime_variation_reference.csv`**. This is the **supported** place to implement the \(r/d\) and x100 chain to match the reference spreadsheet.
- **`calculate_s_prime_from_params`** consumes **already-derived** asymptotes and EC50 only.
- **Loaders** (`SPrime.load`, `RawDataset.load_from_file`, `load_from_dataframe`, etc.): require keyword **`response_normalization=`** (`asymptote_normalized` \| `response_scale`) at import; enforce **`Control_Response`** when **`skip_control_response_normalization=False`** and raw data are present. **`process`** applies test/control ratio (when strict) then the chosen **`response_pipeline`** composable. Keep [development.md](../development.md) in sync when behavior changes.

---

## 9. Glossary (quick reference)

| Phrase | Meaning |
|--------|--------|
| **Raw path** | Concentration + response series -> fit -> parameters -> S'. |
| **Pre-calculated path** | EC50 + asymptotes (etc.) supplied -> S' without refitting. |
| **DMSO / vehicle normalization** | \(r/d\); first sense of "normalized." |
| **Asymptote-normalized / max normalization** | Second sense of "normalized"; scale ratios so **max = 1** (then ×100) before fit. |
| **Response-scale** | No max normalization; historical comparability. |
| **x100** | Conventional scale factor; not biologically special by itself. |
| **`skip_control_response_normalization`** | Opt out of **strict `Control_Response` validation** at import; default **False**. Does **not** (today) turn on automatic \(r/d\) math--that is **`response_pipeline`**. |
| **`Control_Response`** | CSV column for per-profile vehicle (**DMSO**) control readout. |

---

## 10. Open integration questions (for product / API design)

These should be resolved when wiring loaders and public APIs:

1. **Propagation** -- How import-time flags (`skip_control_response_normalization`, future `s_prime_method`) are stored on **`DoseResponseProfile`** / metadata and surfaced in **`process()`** / exports.
2. **Validation rules** -- When **`Control_Response`** is missing: hard fail vs warning; interaction with **`skip_control_response_normalization`** only.
3. **Required `s_prime_method`** -- Tied to (1); whether `asymptote_normalized` vs `response_scale` is **mandatory** on `load()` / `process()` or logged defaults.
4. **Pre-calculated + raw both present** -- Covered by **`allow_overwrite_precalc_params`** (overwrites pre-calculated EC50, asymptotes, Hill_Slope, r^2 when refitting from raw); document interaction with normalization flags once loaders apply preprocessing.
5. **Reporting** -- Record which branches ran for reproducibility.

---

*End of document.*
