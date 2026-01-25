# Basic Usage Guide

This guide walks you through using sprime to analyze high-throughput screening data, from loading raw data to calculating delta S' values for comparative analysis.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Data Format](#data-format)
3. [Loading Data](#loading-data)
4. [Processing Data](#processing-data)
5. [Calculating Delta S'](#calculating-delta-s)
6. [Complete Example](#complete-example)
7. [Exporting Results](#exporting-results)
8. [Data Quality Reporting](#data-quality-reporting)
9. [Advanced Usage](#advanced-usage)
10. [Running the Test Suite](#running-the-test-suite)
11. [Troubleshooting](#troubleshooting)

## Quick Start

The basic workflow in sprime is:

1. **Load** raw data from CSV → `RawDataset`
2. **Process** data (fit curves, calculate S') → `ScreeningDataset`
3. **Analyze** using delta S' for comparative analysis

```python
from sprime import SPrime as sp

# Load and process data (use your CSV, or download a sample from the README Quick Start)
raw_data, _ = sp.load("your_data.csv")
screening_data, _ = sp.process(raw_data)

# Calculate delta S' for comparative analysis
delta_results = screening_data.calculate_delta_s_prime(
    reference_cell_lines="normal_cell_line",
    test_cell_lines=["tumor_cell_line"]
)
```

## Data Format

sprime expects CSV files with the following structure:

### Required Columns

- **Compound Name**: Name of the compound/drug
- **Compound_ID**: Unique identifier for the compound (required)
- **Cell_Line**: Name of the cell line being tested
- **Concentration_Units**: Required when using raw dose-response data (Path A). Units for concentration values (e.g. `microM`, `nM`). See [Supported concentration units](#supported-concentration-units) below.

### Optional Columns

- **pubchem_sid**: PubChem substance identifier
- **SMILES**: Chemical structure notation
- **Cell_Line_Ref_ID**: Reference identifier for the cell line
- **NCGCID**: Optional pass-through per compound (not used for validation)

### Path A: columns vs list

Raw data can be in two layouts. Use **columns** (default) or **list** via `values_as`:

- **`values_as="columns"`** (default): One column per value — `Data0`..`DataN`, `Conc0`..`ConcN`.
- **`values_as="list"`**: Two columns — `Responses` and `Concentrations`. Each cell holds comma-separated values (e.g. `"4000,300,2"`, `"10,3,0.1"`). Order must match; same length; ≥4 pairs. If your CSV is comma-delimited, quote those cells.

```python
raw_data, _ = sp.load("your_data.csv", values_as="columns")  # default
# or
raw_data, _ = sp.load("your_data_list_format.csv", values_as="list")
```

**Templates:** [template_raw.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_raw.csv) (columns), [template_raw_list.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_raw_list.csv) (list).

### Dose-Response Data Columns

For raw dose-response data (**columns** format), use columns named:
- **Data0, Data1, Data2, ... DataN**: Response values at each concentration
- **Conc0, Conc1, Conc2, ... ConcN**: Corresponding concentration values

For **list** format, use **Responses** and **Concentrations** (comma-separated values in one cell each).

### Supported concentration units

When using raw data (Path A), `Concentration_Units` is required. All values are converted to **microM** internally. Supported units (case-insensitive), smallest to largest: **fM** (`fm`, `femtom`); **pM** (`pm`, `picom`); **nM** (`nm`, `nanom`); **microM** (`µM`, `um`, `microm`, `micro`); **mM** (`mm`, `millim`); **M** (`m`, `mol`).

### Pre-calculated Parameters (Optional)

If you already have fitted Hill curve parameters:
- **AC50** or **ec50**: Half-maximal concentration
- **Upper** or **Infinity**: Upper asymptote
- **Lower** or **Zero**: Lower asymptote
- **Hill_Slope** (or **Hill**, **slope**): Hill coefficient
- **r2** or **R²**: R-squared goodness of fit

### Rows are literal

Each row is taken as-is. Empty values are treated as null; there is no forward-filling from previous or subsequent rows. If a cell is empty (e.g. MOA), that field is null for that row. Provide explicit values in every row where you want them.

**Template files:** [template_raw.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_raw.csv) (raw columns), [template_raw_list.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_raw_list.csv) (raw list), [template_precalc.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_precalc.csv) (pre-calculated). Templates put required columns first. Your CSV must use the same header names; column order in the file does not matter.

## Loading Data

### From CSV File

Use your own file path, or [download a sample](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_delta.csv) (raw) / [demo_data_precalc.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_precalc.csv) (pre-calculated), then load:

```python
from sprime import SPrime as sp

# You can download samples from the URLs above, then load your file.
raw_data, _ = sp.load("your_data.csv")

# Or specify an assay name
raw_data, _ = sp.load("your_data.csv", assay_name="HTS001")
```

### Inspecting Loaded Data

```python
# Check how many profiles were loaded
print(f"Loaded {len(raw_data)} profiles")

# Iterate over profiles
for profile in raw_data.profiles:
    print(f"{profile.compound.name} vs {profile.cell_line.name}")
    if profile.concentrations:
        print(f"  Has raw data: {len(profile.concentrations)} data points")
    if profile.hill_params:
        print(f"  Has pre-calculated params: EC50={profile.hill_params.ec50}")
    if profile.metadata:
        print(f"  Metadata: {profile.metadata}")
```

### Getting Specific Profiles

```python
# Get a specific profile by compound and cell line
profile = raw_data.get_profile("DRUG001", "Cell_Line_1")
if profile:
    print(f"Found profile: {profile.compound.name}")
else:
    print("Profile not found")

# You can also pass Compound and CellLine objects
from sprime import Compound, CellLine
compound = Compound(name="Drug A", drug_id="DRUG001")
cell_line = CellLine(name="Cell_Line_1")
profile = raw_data.get_profile(compound, cell_line)
```

### Working with Assay Information

```python
# Access assay information
print(f"Assay name: {raw_data.assay.name}")
print(f"Assay description: {raw_data.assay.description}")

# Create assay with additional metadata
from sprime import Assay
assay = Assay(
    name="HTS002",
    description="High-throughput screen for neurofibroma cell lines",
    screen_id="HTS002",
    readout_type="activity",
    time_profile="48Hr"
)
```

## Processing Data

Processing fits Hill curves (if needed) and calculates S' values:

```python
from sprime import SPrime as sp

# Load data (use your CSV, or download a sample from the README Quick Start)
raw_data, _ = sp.load("your_data.csv")

# Process: fit curves and calculate S'
screening_data, _ = sp.process(raw_data)

# Access results
for profile in screening_data.profiles:
    print(f"{profile.compound.name}: S' = {profile.s_prime:.2f}")
```

### Processing option: `allow_overwrite_hill_coefficients`

When your CSV has **both** raw dose-response data (DATA*/CONC*) **and** pre-calculated Hill params (AC50, Upper, Lower, Hill_Slope, r2), sprime fits from raw and would overwrite those pre-calc values. By default it **raises** unless you explicitly allow overwriting:

```python
# CSV has both raw DATA/CONC and AC50/Upper/Lower columns – allow overwrite
screening_data, _ = sp.process(raw_data, allow_overwrite_hill_coefficients=True)
```

- **`allow_overwrite_hill_coefficients=False`** (default): Raise if we would overwrite pre-calc Hill params.
- **`allow_overwrite_hill_coefficients=True`**: Fit from raw, overwrite pre-calc, and **log a warning** (console + report) that pre-calc was overwritten.

Use `True` when you intentionally want to refit from raw (e.g. demo data, legacy files with both). Same option exists for `get_s_prime_from_data` and `get_s_primes_from_file`.

### Custom Fitting Parameters

You can pass **curve-fitting** parameters (distinct from `allow_overwrite_hill_coefficients`) to control the Hill fit:

```python
# Curve-fitting parameters (maxfev, initial_ec50, etc.)
screening_data, _ = sp.process(
    raw_data,
    curve_direction="up",  # Force increasing curve
    maxfev=10000,          # Faster fitting
    initial_ec50=10.0      # Better initial guess
)

# When CSV has both raw + pre-calc, allow overwrite and optionally pass fit params
screening_data, _ = sp.process(
    raw_data,
    allow_overwrite_hill_coefficients=True,
    maxfev=10000
)
```

**Common fitting parameters:**
- `curve_direction`: `"up"` (increasing), `"down"` (decreasing), or `None` (auto-detect)
- `maxfev`: Maximum function evaluations (default: 3,000,000)
- `initial_lower`, `initial_upper`, `initial_ec50`, `initial_hill_coefficient`: Initial parameter guesses
- `bounds`: Parameter bounds as `([lower_bounds], [upper_bounds])` tuples
- `zero_replacement`: Value to replace zero concentrations (default: 1e-24)

See [Hill Curve Fitting Configuration](hill_curve_fitting_configuration.md) for all available parameters.

### Processing Individual Profiles

You can also process profiles individually:

```python
# Get a profile from raw data
profile = raw_data.get_profile("DRUG001", "Cell_Line_1")

# Fit curve and calculate S' in one step
s_prime = profile.fit_and_calculate_s_prime()

# Or do it step by step
profile.fit_hill_curve()
s_prime = profile.calculate_s_prime()

# Access fitted parameters
print(f"EC50: {profile.hill_params.ec50}")
print(f"Upper: {profile.hill_params.upper}")
print(f"Lower: {profile.hill_params.lower}")
print(f"R²: {profile.hill_params.r_squared}")
print(f"S': {profile.s_prime}")
```

### Handling Processing Errors

If processing fails for some profiles, you can handle errors gracefully:

```python
from sprime import SPrime as sp

raw_data, _ = sp.load("your_data.csv")
screening_data = ScreeningDataset(assay=raw_data.assay)

# Process profiles one at a time to handle errors
for profile in raw_data.profiles:
    try:
        # Create a copy for processing
        processed = DoseResponseProfile(
            compound=profile.compound,
            cell_line=profile.cell_line,
            assay=profile.assay,
            concentrations=profile.concentrations,
            responses=profile.responses,
            hill_params=profile.hill_params,
            s_prime=profile.s_prime
        )
        
        # Process if needed
        if processed.hill_params is None and processed.concentrations:
            processed.fit_hill_curve()
        if processed.s_prime is None:
            processed.calculate_s_prime()
        
        screening_data.add_profile(processed)
    except (ValueError, RuntimeError) as e:
        print(f"Failed to process {profile.compound.name} vs {profile.cell_line.name}: {e}")
        continue
```

## Calculating Delta S'

Delta S' (ΔS') compares drug responses between reference and test cell lines. This is useful for identifying compounds with selective activity.

### Basic Delta S' Calculation

```python
# Calculate delta S' = S'(reference) - S'(test)
delta_results = screening_data.calculate_delta_s_prime(
    reference_cell_lines="normal_tissue",
    test_cell_lines=["tumor_cell_line"]
)

# Access results
for ref_cellline, comparisons in delta_results.items():
    print(f"\nReference: {ref_cellline}")
    for comp in comparisons:
        print(f"  {comp['compound_name']}: Delta S' = {comp['delta_s_prime']:.2f}")
```

### Multiple Reference and Test Cell Lines

```python
# Compare multiple cell lines
delta_results = screening_data.calculate_delta_s_prime(
    reference_cell_lines=["normal_tissue_1", "normal_tissue_2"],
    test_cell_lines=["tumor_line_1", "tumor_line_2", "tumor_line_3"]
)
```

### Understanding Delta S' Results

Delta S' results are organized by reference cell line:

```python
delta_results = {
    "reference_cell_line_1": [
        {
            "compound_name": "Drug A",
            "drug_id": "DRUG001",
            "test_cell_line": "tumor_line_1",
            "s_prime_reference": 2.5,
            "s_prime_test": 4.0,
            "delta_s_prime": -1.5  # Negative = more effective in test
        },
        # ... more comparisons
    ],
    "reference_cell_line_2": [
        # ... comparisons for second reference
    ]
}
```

**Interpreting Delta S':**
- **Negative values**: More effective in test cell line (higher S' in test)
- **Positive values**: More effective in reference cell line (higher S' in reference)
- **More negative = more selective**: Compounds with more negative delta S' are more selective for the test (e.g., tumor) cell line

### Ranking by Delta S'

```python
# Calculate delta S'
delta_results = screening_data.calculate_delta_s_prime(
    reference_cell_lines="normal_tissue",
    test_cell_lines=["tumor_line"]
)

# Extract and sort by delta S' (most negative = most selective)
for ref_cellline, comparisons in delta_results.items():
    sorted_comps = sorted(comparisons, key=lambda x: x['delta_s_prime'])
    
    print(f"\nRanking for {ref_cellline} vs tumor_line:")
    for rank, comp in enumerate(sorted_comps, start=1):
        print(f"{rank}. {comp['compound_name']}: Delta S' = {comp['delta_s_prime']:.2f}")
```

## Complete Example

Here's a complete example using the demo data:

```python
from sprime import SPrime as sp

# 1. Load data (use your CSV, or download docs/usage/demo_data_s_prime.csv from the repo)
print("Loading data...")
raw_data, _ = sp.load("docs/usage/demo_data_s_prime.csv")
print(f"Loaded {len(raw_data)} profiles")

# 2. Process data (fit curves, calculate S')
print("\nProcessing data (fitting curves, calculating S')...")
screening_data, _ = sp.process(raw_data)
print("Processing complete!")

# 3. Display S' values
print("\n=== S' Values ===")
for profile in screening_data.profiles:
    print(f"{profile.compound.name} vs {profile.cell_line.name}: "
          f"S' = {profile.s_prime:.2f}")

# 4. Calculate delta S' for comparative analysis
print("\n=== Delta S' Analysis ===")
delta_results = screening_data.calculate_delta_s_prime(
    reference_cell_lines="Normal_Tissue",
    test_cell_lines=["Tumor_Cell_Line"]
)

# 5. Display and rank results
for ref_cellline, comparisons in delta_results.items():
    print(f"\nReference: {ref_cellline}")
    
    # Sort by delta S' (most negative = most selective for tumor)
    sorted_comps = sorted(comparisons, key=lambda x: x['delta_s_prime'])
    
    print(f"\nRanking (most selective for tumor first):")
    for rank, comp in enumerate(sorted_comps, start=1):
        print(f"  {rank}. {comp['compound_name']}: "
              f"Delta S' = {comp['delta_s_prime']:.2f} "
              f"(S' ref={comp['s_prime_reference']:.2f}, "
              f"S' test={comp['s_prime_test']:.2f})")
```

## Exporting Results

### Export to CSV

The easiest way to export results is using the built-in CSV export methods:

```python
# Export all profiles to CSV
screening_data.export_to_csv("master_s_prime_table.csv", include_metadata=True)

# Export delta S' results to CSV
delta_results = screening_data.calculate_delta_s_prime(...)
ScreeningDataset.export_delta_s_prime_to_csv(delta_results, "delta_s_prime_table.csv")
```

The `export_to_csv()` method includes all profile information including Hill curve parameters, S' values, ranking, and optional metadata (any non-reserved columns from your CSV). Delta S' export adds reserved compound-level columns MOA and drug targets (resolved from common header variants) plus any optional headings you specify.

### Export to List of Dictionaries

For programmatic access or custom export formats:

```python
# Export all profiles as dictionaries
results = screening_data.to_dict_list()

# Save to JSON, CSV, etc.
import json
with open("results.json", "w") as f:
    json.dump(results, f, indent=2)
```

### Working with In-Memory Data

If you have data as a list of dictionaries (e.g., from a database query), you can process it directly:

```python
from sprime import get_s_prime_from_data, calculate_delta_s_prime

# List of dicts matching CSV row format
list_of_rows = [
    {
        'Compound Name': 'Drug A',
        'Compound_ID': 'DRUG001',
        'Cell_Line': 'Cell_Line_1',
        'Data0': '10', 'Data1': '20', 'Data2': '50', 'Data3': '90',
        'Conc0': '0.1', 'Conc1': '1', 'Conc2': '10', 'Conc3': '100',
    },
    # ... more rows
]

# Calculate S' values directly from list of dicts
results = get_s_prime_from_data(list_of_rows)

# If rows have both raw DATA/CONC and pre-calc AC50/Upper/Lower, allow overwrite:
# results = get_s_prime_from_data(list_of_rows, allow_overwrite_hill_coefficients=True)

# Calculate delta S' from list of dicts
delta_results = calculate_delta_s_prime(
    results,  # Can use list of dicts directly
    reference_cell_line_names="Reference_Cell_Line",
    test_cell_line_names="Test_Cell_Line"
)
```

## Next Steps

- Learn about [Hill Curve Fitting Configuration](hill_curve_fitting_configuration.md)
- Read about [Understanding 4PL Dose-Response Curves](../README_4PL_Dose_Response.md)
- Review [Background and Concepts](../background_and_concepts.md) for qHTS terminology and S' metric details

## Troubleshooting

### "Pre-calculated Hill parameters would be overwritten" Error

Your CSV has **both** raw dose-response columns (Data0..DataN, Conc0..ConcN) **and** pre-calc Hill params (AC50, Upper, Lower, Hill_Slope, r2). By default sprime raises to avoid silently overwriting user-supplied values.

**Fix:** Set `allow_overwrite_hill_coefficients=True` when you intend to refit from raw and overwrite pre-calc:

```python
screening_data, _ = sp.process(raw_data, allow_overwrite_hill_coefficients=True)
# or
results = get_s_prime_from_data(list_of_rows, allow_overwrite_hill_coefficients=True)
```

When you allow overwrite, sprime logs a **warning** (console + report) each time pre-calc Hill params are overwritten by fitted values.

### "Cannot fit curve" Error

This means a profile has no raw data and no pre-calculated Hill parameters. Ensure your CSV has either:
- Data0..DataN and Conc0..ConcN columns (raw data), OR
- AC50, Upper, Lower columns (pre-calculated parameters)

**Common causes:**
- Missing concentration or response columns
- Empty or malformed numeric values
- Insufficient data points (need at least 4 points for 4-parameter fit)

### Low R² Values

If curve fits have low R² (< 0.7), consider:
- Checking data quality (outliers, missing points)
- Ensuring concentrations span the full response range
- Adjusting initial parameter guesses
- See [Hill Curve Fitting Configuration](hill_curve_fitting_configuration.md) for options

**Note:** Low R² values may indicate:
- Non-sigmoidal dose-response relationship
- Experimental artifacts or errors
- Insufficient concentration range

### No Delta S' Results

If `calculate_delta_s_prime()` returns empty results:
- Verify cell line names match exactly (case-sensitive)
- Ensure you have data for both reference and test cell lines
- Check that S' values were calculated successfully

**Debugging tip:** Check individual profiles first:
```python
profile = screening_data.get_profile("DRUG001", "Cell_Line_1")
if profile:
    print(f"S' = {profile.s_prime}")
else:
    print("Profile not found")
```

### Import Errors

If you see `ImportError: Hill curve fitting requires scipy`:
- Install scipy: `pip install scipy`
- Note: numpy is automatically installed as a dependency of scipy

### CSV Export Issues

If CSV export fails:
- Ensure the output directory exists
- Check file permissions
- Verify the ScreeningDataset has profiles with S' values calculated

### Metadata Not Appearing

If metadata is not appearing in exports:
- Ensure non-reserved columns exist in your CSV (e.g. MOA, drug targets, Target, MoA); these are stored as generic metadata under the exact header.
- Check that `include_metadata=True` is set when calling `export_to_csv()`.
- For delta S' tables, MOA and drug targets are reserved compound-level columns and are resolved from common variants (MOA/MoA/moa, drug targets/Target/target) at export time.
- Verify metadata was extracted during loading (check `profile.metadata`).

## Advanced Usage

### Working with Multiple Assays

If you have data from multiple assays, process them separately:

```python
from sprime import SPrime as sp

# Load and process each assay separately
assay1_data, _ = sp.load("assay1_data.csv", assay_name="HTS001")
assay1_processed, _ = sp.process(assay1_data)

assay2_data, _ = sp.load("assay2_data.csv", assay_name="HTS002")
assay2_processed, _ = sp.process(assay2_data)

# Compare results across assays (be cautious - different conditions)
for profile1 in assay1_processed.profiles:
    profile2 = assay2_processed.get_profile(
        profile1.compound.drug_id,
        profile1.cell_line.name
    )
    if profile2:
        print(f"{profile1.compound.name}: S' in HTS001={profile1.s_prime:.2f}, "
              f"S' in HTS002={profile2.s_prime:.2f}")
```

### Filtering and Querying Results

```python
# Filter profiles by S' threshold
high_s_prime = [
    p for p in screening_data.profiles
    if p.s_prime and p.s_prime > 2.0
]

# Filter by compound
from sprime import Compound
target_compound = Compound(name="Target Drug", drug_id="TARGET001")
target_profiles = [
    p for p in screening_data.profiles
    if p.compound.drug_id == target_compound.drug_id
]

# Filter by cell line
tumor_profiles = [
    p for p in screening_data.profiles
    if "tumor" in p.cell_line.name.lower()
]

# Sort by S' value
sorted_profiles = sorted(
    screening_data.profiles,
    key=lambda p: p.s_prime if p.s_prime else float('-inf'),
    reverse=True
)
```

### Creating dose-response profiles from scratch

You can create `DoseResponseProfile` objects programmatically (without loading from CSV), then fit and calculate S':

```python
from sprime import DoseResponseProfile, Compound, CellLine, Assay

compound = Compound(name="Trifluoperazine", drug_id="NCGC00013226-15")
cell_line = CellLine(name="ipNF96.11C")
assay = Assay(name="HTS002", readout_type="activity")

profile = DoseResponseProfile(
    compound=compound,
    cell_line=cell_line,
    assay=assay,
    concentrations=[1.30e-9, 3.91e-9, 1.17e-8, 3.52e-8, 1.06e-7, 3.17e-7, 9.50e-7, 2.85e-6, 8.55e-6, 2.56e-5, 7.69e-5],
    responses=[-63.23, -66.40, -67.04, -65.47, -63.59, -60.30, -47.43, 46.75, 69.12, 97.97, 85.27],
    concentration_units="microM"
)

s_prime = profile.fit_and_calculate_s_prime()
print(f"S' = {s_prime:.2f}")
```

### Working with Pre-calculated Hill Parameters

If your data already has fitted parameters, you can skip curve fitting:

```python
from sprime import SPrime as sp, HillCurveParams

# Load data with pre-calculated parameters
raw_data, _ = sp.load("data_with_params.csv")

# Profiles will have hill_params set
for profile in raw_data.profiles:
    if profile.hill_params:
        # Just calculate S' - no fitting needed
        s_prime = profile.calculate_s_prime()
        print(f"{profile.compound.name}: S' = {s_prime:.2f}")

# Or manually set parameters
profile = raw_data.get_profile("DRUG001", "Cell_Line_1")
if profile:
    profile.hill_params = HillCurveParams(
        ec50=10.0,
        upper=100.0,
        lower=0.0,
        hill_coefficient=1.5,
        r_squared=0.95
    )
    s_prime = profile.calculate_s_prime()
```

### Batch Processing Multiple Files

```python
from pathlib import Path
from sprime import SPrime as sp

# Process all CSV files in a directory
data_dir = Path("screening_data")
results = []

for csv_file in data_dir.glob("*.csv"):
    print(f"Processing {csv_file.name}...")
    try:
        raw_data, _ = sp.load(csv_file)
        screening_data, _ = sp.process(raw_data)
        results.append({
            'file': csv_file.name,
            'assay': screening_data.assay.name,
            'num_profiles': len(screening_data),
            'data': screening_data
        })
    except Exception as e:
        print(f"Error processing {csv_file.name}: {e}")
        continue

# Combine results or analyze separately
for result in results:
    print(f"{result['file']}: {result['num_profiles']} profiles")
```

## Data Quality Reporting

sprime automatically tracks data quality issues and warnings during data loading and processing. By default, a summary report is printed to the console when processing data. You can configure reporting to disable console output, enable verbose mode, or write detailed log files.

### Default Behavior

When you load or process data, sprime automatically prints a summary report to the console:

```python
from sprime import SPrime as sp

# Console summary is printed automatically
raw_data, load_report = sp.load("data.csv")
screening_data, process_report = sp.process(raw_data)
```

**Example Console Output:**
```
============================================================
DATA PROCESSING SUMMARY
============================================================
Total Rows:                 247
Rows Processed:             245
Rows Skipped:               2
Compounds Loaded:           58
Profiles Created:           233
Profiles with S' Calculated: 221
Profiles Failed:            12

DATA QUALITY ISSUES:
  Missing Drug IDs:          2
  Missing Compound Names:    5
  Insufficient Data Points:  8
  Invalid Numeric Values:    3

WARNINGS: 25 total
  MISSING_DATA: 10
    Row 12: Missing Compound_ID, row skipped (CL:A549)
    Row 23: Insufficient data points: 3 found (ID:DRUG004, CL:HeLa)
    ... and 8 more (see log file for details)
============================================================
```

### Configuring Reporting

Use `ReportingConfig` to control reporting behavior globally:

```python
from sprime import SPrime as sp, ReportingConfig, ConsoleOutput

# Configure reporting settings
ReportingConfig.configure(
    log_to_file=True,                    # Enable log file writing
    log_filepath="processing.log",       # Optional: custom log file path
    console_output=ConsoleOutput.SUMMARY  # Console verbosity
)
```

### Console Output Options

**Summary Mode (Default):**
```python
from sprime import ReportingConfig, ConsoleOutput

# Default: Brief summary with counts and first 3 examples per category
ReportingConfig.configure(console_output=ConsoleOutput.SUMMARY)
```

**Verbose Mode:**
```python
# Print all warnings to console with full details
ReportingConfig.configure(console_output=ConsoleOutput.VERBOSE)
```

**Disable Console Output:**
```python
# No console output (useful for batch processing)
ReportingConfig.configure(console_output=ConsoleOutput.NONE)
```

### Log File Writing

Enable detailed log file writing:

```python
from sprime import SPrime as sp, ReportingConfig

# Enable log file (auto-generated from input filename)
ReportingConfig.configure(log_to_file=True)

raw_data, report = sp.load("data.csv")
# → Creates "data_processing.log" automatically
```

**Custom Log File Path:**
```python
from sprime import SPrime as sp, ReportingConfig

# Specify custom log file path
ReportingConfig.configure(
    log_to_file=True,
    log_filepath="my_custom_log.log"
)

raw_data, report = sp.load("data.csv")
# → Writes to "my_custom_log.log"
```

**Log File Format:**
The log file contains:
- Summary metrics (rows processed, compounds loaded, profiles created)
- Data quality issue counts
- Detailed warnings grouped by category with row numbers and context
- Each warning includes: Drug ID, Compound Name, Cell Line, Field Name

### Complete Examples

**Example 1: Default Usage (Console Summary Only)**
```python
from sprime import SPrime as sp

# Default: Console summary printed, no log file
raw_data, load_report = sp.load("data.csv")
screening_data, process_report = sp.process(raw_data)
```

**Example 2: Enable Log File**
```python
from sprime import SPrime as sp, ReportingConfig

# Enable log file writing
ReportingConfig.configure(log_to_file=True)

raw_data, load_report = sp.load("data.csv")
# → Console summary printed
# → Log file written: "data_processing.log"

screening_data, process_report = sp.process(raw_data)
# → Console summary printed
# → Log file appended/updated
```

**Example 3: Verbose Console Output**
```python
from sprime import SPrime as sp, ReportingConfig, ConsoleOutput

# Verbose console output
ReportingConfig.configure(console_output=ConsoleOutput.VERBOSE)

raw_data, report = sp.load("data.csv")
# → All warnings printed to console with full details
```

**Example 4: Silent Mode (Log File Only)**
```python
from sprime import SPrime as sp, ReportingConfig, ConsoleOutput

# No console output, but write log file
ReportingConfig.configure(
    log_to_file=True,
    console_output=ConsoleOutput.NONE
)

raw_data, report = sp.load("data.csv")
# → No console output
# → Log file written: "data_processing.log"
```

**Example 5: Reset to Defaults**
```python
from sprime import ReportingConfig

# Reset all settings to defaults
ReportingConfig.reset()
# → Console summary enabled, log file disabled
```

### Understanding Warnings

**Warning Categories:**
- `MISSING_DATA`: Missing required fields (Compound_ID, Cell_Line, insufficient data points)
- `DATA_QUALITY`: Invalid values, non-numeric data, missing optional fields
- `NUMERICAL`: NaN/Inf values encountered
- `CURVE_FIT`: Fitting failures, poor fit quality (R² < 0.7)
- `CALCULATION`: S' calculation failures

**Row Numbers:**
- Row numbers start at 2 (row 1 is the CSV header)
- Row number 0 indicates warnings from processing (not from CSV loading)
- Fully blank rows are skipped silently (not logged)

**Using Reports Programmatically:**
```python
from sprime import SPrime as sp

raw_data, report = sp.load("data.csv")

# Access report metrics
print(f"Compounds loaded: {report.compounds_loaded}")
print(f"Warnings: {len(report.warnings)}")

# Access individual warnings
for warning in report.warnings:
    if warning.category == "MISSING_DATA":
        print(f"Row {warning.row_number}: {warning.message}")
```

### Best Practices

1. **Development/Interactive Use**: Use default settings (console summary) to see issues immediately
2. **Batch Processing**: Disable console output (`ConsoleOutput.NONE`) and enable log files
3. **Debugging**: Use verbose mode (`ConsoleOutput.VERBOSE`) to see all warnings in console
4. **Production**: Enable log files for audit trails and quality tracking

## Running the Test Suite

sprime includes a comprehensive test suite to verify functionality and catch regressions. Running the tests ensures that the library is working correctly in your environment.

### Prerequisites

Install the development dependencies:

```bash
pip install -e ".[dev]"
```

Or install pytest directly:

```bash
pip install pytest pytest-cov
```

### Running All Tests

Run the complete test suite:

```bash
# From the project root directory
pytest tests/
```

This will run all tests in the `tests/` directory:
- `test_hill_fitting.py` - Tests for Hill curve fitting functionality
- `test_sprime.py` - Tests for core sprime module functionality
- `test_integration.py` - End-to-end integration tests

### Running Specific Test Files

Run tests for a specific module:

```bash
# Test only Hill curve fitting
pytest tests/test_hill_fitting.py

# Test only core sprime functionality
pytest tests/test_sprime.py

# Test only integration workflows
pytest tests/test_integration.py
```

### Running Specific Test Classes or Functions

Run a specific test class or function:

```bash
# Run a specific test class
pytest tests/test_sprime.py::TestDoseResponseProfile

# Run a specific test function
pytest tests/test_sprime.py::TestDoseResponseProfile::test_fit_hill_curve

# Run tests matching a pattern
pytest tests/ -k "delta"  # Runs all tests with "delta" in the name
```

### Verbose Output

Get more detailed output:

```bash
# Verbose output showing each test
pytest tests/ -v

# Very verbose output with print statements
pytest tests/ -v -s

# Show local variables on failure
pytest tests/ -v -l
```

### Test Coverage

Generate a coverage report to see which code is tested:

```bash
# Run tests with coverage
pytest tests/ --cov=src/sprime --cov-report=html

# View coverage report
# Open htmlcov/index.html in your browser

# Or get a terminal report
pytest tests/ --cov=src/sprime --cov-report=term
```

### Common Test Scenarios

The test suite covers:

**Core Functionality:**
- Value object creation (Compound, CellLine, Assay, HillCurveParams)
- Dose-response profile operations
- RawDataset loading from CSV
- ScreeningDataset processing and analysis
- S' value calculations
- Delta S' calculations

**Edge Cases:**
- Missing or malformed data
- Insufficient data points
- Failed curve fits
- Missing profiles in delta S' calculations
- CSV parsing with various formats

**Integration Tests:**
- Complete workflows (Load → Process → Analyze)
- Multiple compounds and cell lines
- Metadata extraction and preservation
- CSV export functionality
- In-memory data processing

### Troubleshooting Tests

**If tests fail:**

1. **Check dependencies**: Ensure scipy and numpy are installed
   ```bash
   pip install scipy numpy
   ```

2. **Check Python version**: sprime requires Python 3.8+
   ```bash
   python --version
   ```

3. **Run with verbose output** to see detailed error messages:
   ```bash
   pytest tests/ -v -s
   ```

4. **Check for import errors**: Ensure you're running from the project root
   ```bash
   cd /path/to/sprime
   pytest tests/
   ```

**Common test failures:**

- **ImportError**: Make sure you're in the project root directory
- **FileNotFoundError**: Some tests create temporary files - ensure write permissions
- **AssertionError**: Check that your data matches expected formats

### Continuous Integration

If contributing to sprime, ensure all tests pass before submitting:

```bash
# Run full test suite
pytest tests/

# Check code style (if configured)
# flake8 src/ tests/

# Run with coverage
pytest tests/ --cov=src/sprime --cov-report=term-missing
```

For more information on testing, see the [pytest documentation](https://docs.pytest.org/).

