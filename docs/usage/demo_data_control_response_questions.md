# `Control_Response` (vehicle / DMSO) -- demo CSV conventions

The templates **`template_raw.csv`** and **`template_raw_list.csv`** include **`Control_Response`** using **`-67`**, matching the **lowest** value in their toy response series (same scale as `DATA*` / `Responses`).

## Two demo families

| Files | Meaning | `Control_Response` | Load with |
|-------|---------|-------------------|-----------|
| **`demo_precontrol_normalized_*`** | Responses already normalized to vehicle / analysis scale (e.g. ipNF-style). | Column **present**, cells **empty** | **`skip_control_response_normalization=True`** (required) |
| **`demo_raw_vehicle_control_*`** | Raw readouts built from the committed reference **`tests/fixtures/SPrime_variation_reference.csv`**: for each 17-AAG row, **final compound concentration uM** -> `CONC*`, **percent nucleus response, mean of triplicate** -> `DATA*`; **`Control_Response` = 35.3** from the **DMSO** row. | **35.3** (per row) | **Default** (`skip_control_response_normalization=False`) |

The loader **validates** `Control_Response` by **default** when the row has a raw curve: value must be non-empty, numeric, and non-zero. Set **`skip_control_response_normalization=True`** only when responses are already control-normalized. The loader does **not** yet apply `sprime.response_pipeline` scaling automatically--use **`pipeline_response_scale`** (or **`pipeline_asymptote_normalized`**) if you need the x100 scale from raw + DMSO.

---

## `demo_precontrol_normalized_s_prime.csv`

- **What it is:** One raw curve (`ipNF96.11C`) plus pre-calculated Hill columns on the same row; **`Control_Response`** is blank.

---

## `demo_precontrol_normalized_raw_list.csv`

- **What it is:** Same biological example as `demo_precontrol_normalized_s_prime.csv`, **list** format (`Responses` / `Concentrations`), one row; **`Control_Response`** blank.

---

## `demo_precontrol_normalized_delta.csv`

- **What it is:** Two rows--same compound, two cell lines (`ipNF96.11C` vs `ipnf05.5 mc`), pre-control-normalized raw curves; **`Control_Response`** blank on both rows.

---

## `demo_precontrol_normalized_precalc.csv`

- **What it is:** Path B only--`AC50`, `Upper`, `Lower`, etc. **`Control_Response`** column is present and blank (uniform schema; ignored for validation on pure pre-calc rows).

---

## `demo_raw_vehicle_control_s_prime.csv` / `_raw_list.csv` / `_delta.csv`

- **What it is:** 17-AAG / NF2-/- / 48h slice from **`SPrime_variation_reference.csv`** (see table above for column mapping). Not hand-made toy numbers--these are the reference sheet values.
- **Delta file:** Same concentrations and **Control_Response**; **Test** line uses responses x0.98 vs **Ref** (synthetic perturbation only for a two-line delta demo--not from the CSV).

---

## After you customize

- For **pre-control-normalized** exports from your pipeline: keep **`Control_Response`** in the header with empty cells and document **`skip_control_response_normalization=True`**.
- For **raw + vehicle**: fill **`Control_Response`** per row and use **default** loader settings (strict), then preprocess responses if your fitting scale requires it.

Update **[Basic usage guide](basic_usage_guide.md)** / **[README](../../README.md)** when you add new demo variants.
