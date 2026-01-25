# sprime

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Development Status](https://img.shields.io/badge/status-alpha-orange)](https://github.com/MoCoMakers/sprime)

**sprime** (S' or S prime) is a Python library for analyzing quantitative high-throughput screening (qHTS) data in preclinical drug discovery studies. The library provides tools for processing dose-response curves, fitting Hill equations, and calculating S' values—a single metric that summarizes a drug's dose-response profile relative to a cell line and assay.

## Overview

sprime enables researchers to:

- **Process qHTS data**: Load and analyze dose-response data from screening campaigns
- **Fit Hill curves**: Automatically fit four-parameter logistic (4PL) models to dose-response data
- **Calculate S' values**: Compute S' (S prime) metrics that combine potency and efficacy into a single score
- **Compare cell lines**: Use delta S' (ΔS') to identify compounds with selective activity between reference and test cell lines
- **Rank compounds**: Systematically prioritize drug candidates based on their dose-response profiles

## About S' (S Prime)

S' is a single value score that summarizes a drug's dose-response curve. The metric is calculated from Hill curve parameters using the formula:

**S' = asinh((Upper - Lower) / EC50)**

This is equivalent to:

$$S' = \ln\left(\frac{\text{Upper} - \text{Lower}}{\text{EC50}} + \sqrt{\left(\frac{\text{Upper} - \text{Lower}}{\text{EC50}}\right)^2 + 1}\right)$$

Where:
- **Upper** = Upper asymptote of the dose-response curve
- **Lower** = Lower asymptote of the dose-response curve  
- **EC50** = Half-maximal effect concentration

Higher S' values indicate stronger drug responses (higher efficacy and/or potency). See [Background and Concepts](docs/background_and_concepts.md#the-s-s-prime-metric) for detailed information.

### Delta S' and Comparative Analysis

**Delta S'** (ΔS') enables quantitative comparison of drug responses between different cell lines within a single assay:

**ΔS' = S'(reference cell line) - S'(test cell line)**

This metric allows researchers to:
- Compare drug effects across cell lines
- Rank compounds by selectivity
- Prioritize drug candidates for further investigation

For detailed information and examples, see [Delta S' for Comparative Analysis](docs/background_and_concepts.md#delta-s-for-comparative-analysis).

This library implements **Generation 2** of the S' methodology, which evolved from the original S metric. See the [Citation](#citation) section below for references.

## Installation

### Requirements

- Python 3.8 or higher
- numpy >= 1.20.0, < 3.0
- scipy >= 1.7.0

### Install from PyPI

```bash
pip install sprime
```

### Install from Source

For development or to get the latest version:

```bash
# Clone the repository
git clone https://github.com/MoCoMakers/sprime.git
cd sprime

# Install in editable mode
pip install -e .

# Or install with development dependencies
pip install -e ".[dev]"
```

### Development Setup

To set up a development environment:

```bash
git clone https://github.com/MoCoMakers/sprime.git
cd sprime
python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -e ".[dev]"
pytest tests/
pre-commit install        # optional: rebuild + stage pdoc_html on commit when .py changes
```

See **[Development guide](docs/development.md)** for full setup, pre-commit, API docs, and CI.

## Quick Start

The basic workflow in sprime is: **Load** raw data from CSV → **Process** (fit curves, calculate S') → **Analyze** (e.g. delta S' for comparative analysis).

**Three ways to load data:**

| Path | `values_as` | Your data | Required columns | What sprime does |
|------|-------------|-----------|------------------|------------------|
| **A — Raw (columns)** | `"columns"` | `DATA0`..`DATAN`, `CONC0`..`CONCN` (one column per value) | `Cell_Line`, `Compound_ID`, `Concentration_Units` | Fits Hill equation coefficients and stats → AC50, Upper, Lower, Hill Slope, r² → S' |
| **A — Raw (list)** | `"list"` | `Responses`, `Concentrations` (comma-separated in one cell each) | `Cell_Line`, `Compound_ID`, `Concentration_Units` | Same as above |
| **B — Pre-calculated** | N/A | `AC50`, `Upper`, `Lower` (optional: Hill Slope, r2, S', etc.) | `Cell_Line`, `Compound_ID` | Uses params as-is; always computes S' (warns if S' was provided) |

Use **Path A (columns)** by default (`values_as="columns"`). Use **Path A (list)** with `sp.load(..., values_as="list")` when your data has `Responses` and `Concentrations` as comma-separated values. For list format, if the CSV is comma-delimited, quote those cells (e.g. `"4000,300,2"`). Any non-reserved columns in your CSV (e.g. MOA, assay_timescale) propagate forward by default into S' output and the master CSV export.

### Path A — Raw data

**Templates:** [template_raw.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_raw.csv) (columns), [template_raw_list.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_raw_list.csv) (list). See [Required headings](#required-headings) below.

The **load** call changes with format; **process** does not. You must pass `values_as="columns"` (default) or `values_as="list"` — format is not auto-detected from headings.

Run from the directory containing the demo CSVs (e.g. `docs/usage/`), or use paths like `docs/usage/demo_data_s_prime.csv`.

**Option 1 — Columns** (DATA0..N, CONC0..N)

```python
from sprime import SPrime as sp

# Download (columns), then save locally:
# https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_s_prime.csv
# https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_delta.csv
# Demo files include both raw DATA/CONC and pre-calc params; refit from raw.

raw_data, _ = sp.load("demo_data_s_prime.csv", values_as="columns")
screening_data, _ = sp.process(raw_data, allow_overwrite_hill_coefficients=True)
screening_data.export_to_csv("master_s_prime_table.csv")
results = screening_data.to_dict_list()
for profile in results:
    print(f"{profile['compound_name']} vs {profile['cell_line']}: S' = {profile['s_prime']:.2f}")
```

**Option 2 — List** (Responses, Concentrations as comma-separated per cell)

```python
from sprime import SPrime as sp

# Download (list), then save locally:
# https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_raw_list.csv

raw_data, _ = sp.load("demo_data_raw_list.csv", values_as="list")
screening_data, _ = sp.process(raw_data)
screening_data.export_to_csv("master_s_prime_table.csv")
results = screening_data.to_dict_list()
for profile in results:
    print(f"{profile['compound_name']} vs {profile['cell_line']}: S' = {profile['s_prime']:.2f}")
```

### Path B — Pre-calculated

**Sample file:** [demo_data_precalc.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_precalc.csv)  
**Template:** [template_precalc.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/template_precalc.csv)

```python
from sprime import SPrime as sp

# Download, then save locally:
# https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_precalc.csv

raw_data, _ = sp.load("demo_data_precalc.csv")
screening_data, _ = sp.process(raw_data)
screening_data.export_to_csv("master_s_prime_table.csv")
results = screening_data.to_dict_list()
for profile in results:
    print(f"{profile['compound_name']} vs {profile['cell_line']}: S' = {profile['s_prime']:.2f}")
```

### Delta S' example

Compare drug responses between reference and test cell lines (e.g. non-tumor vs tumor). Delta S' = S'(reference) − S'(test); more negative = more effective in test tissue.

**Demo CSV (two cell lines):** [demo_data_delta.csv](https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_delta.csv)

```python
from sprime import SPrime as sp, ScreeningDataset

# Download, then save locally:
# https://raw.githubusercontent.com/MoCoMakers/sprime/refs/heads/main/docs/usage/demo_data_delta.csv

raw_data, _ = sp.load("demo_data_delta.csv", values_as="columns")
screening_data, _ = sp.process(raw_data, allow_overwrite_hill_coefficients=True)
screening_data.export_to_csv("master_s_prime_table.csv")

delta_results = screening_data.calculate_delta_s_prime(
    reference_cell_lines=["ipnf05.5 mc"],   # e.g. non-tumor; list supports multiple
    test_cell_lines=["ipNF96.11C"]          # e.g. tumor; list supports multiple
)
# Opt-in: add compound-level columns (1:1 in ref & test) to delta output and export:
#   headings_one_to_one_in_ref_and_test=["assay_timescale"],  # etc.
#   source_profile="test",   # "ref" or "test" — which profile to use for those values
# Pass the same list to export_delta_s_prime_to_csv(..., headings_one_to_one_in_ref_and_test=[...])
ScreeningDataset.export_delta_s_prime_to_csv(delta_results, "delta_s_prime_table.csv")

for ref_cl, comparisons in delta_results.items():
    for c in comparisons:
        print(f"{c['compound_name']}: Delta S' = {c['delta_s_prime']:.2f}")
```

### Required headings

- **All paths:** `Cell_Line`; and `Compound_ID`. (`NCGCID` is optional pass-through per compound.)
- **Path A (raw):** **`Concentration_Units`** (required). Either **(columns)** `DATA0`..`DATA N`, `CONC0`..`CONC N` (same N), or **(list)** `Responses` and `Concentrations` (comma-separated values in one cell each). Use `values_as="columns"` (default) or `values_as="list"`. See [Supported concentration units](#supported-concentration-units) below.
- **Path B (pre-calc):** `AC50` (or `ec50`), `Upper` (or `Infinity`), `Lower` (or `Zero`). Optional: `Hill_Slope`, `r2`, `S'`.

Template files list the exact headers (required columns first); your CSV should match those. Column order in the file does not matter—matching is by header name only.

### Supported concentration units

For Path A (raw), `Concentration_Units` must be present. All values are converted to **microM** internally. Supported units (case-insensitive), smallest to largest: **fM** (`fm`, `femtom`); **pM** (`pm`, `picom`); **nM** (`nm`, `nanom`); **microM** (`µM`, `um`, `microm`, `micro`); **mM** (`mm`, `millim`); **M** (`m`, `mol`).

### Next steps

- **Format details:** [Basic Usage Guide](docs/usage/basic_usage_guide.md)
- **Run all scenarios:** [Demo](docs/usage/demo.py)
- **Individual profiles & pre-calculated parameters:** [Basic Usage Guide](docs/usage/basic_usage_guide.md) (*Processing Individual Profiles*, *Creating dose-response profiles from scratch*, *Working with Pre-calculated Hill Parameters*)

## Key Features

- **Automatic Curve Fitting**: Fits four-parameter logistic (Hill) curves to dose-response data
- **Flexible Input**: Supports raw dose-response data or pre-calculated Hill parameters. When both exist, use `allow_overwrite_hill_coefficients=True` to refit from raw; overwrites are logged as warnings.
- **CSV Loading**: Handles common screening data formats; rows are literal (empty = null, no forward-filling)
- **Delta S' Analysis**: Compare drug responses across cell lines within a single assay
- **CSV Export**: Built-in methods to export results and delta S' comparisons to CSV
- **In-Memory Processing**: Process data directly from list of dictionaries without CSV files
- **Metadata Support**: Extracts and preserves metadata (MOA, drug targets) from CSV files
- **Type-Safe**: Uses Python dataclasses and type hints throughout
- **Comprehensive Documentation**: Includes guides for users and detailed technical documentation

## Documentation

### Getting Started

- **[API Reference](https://mocomakers.github.io/sprime/)** - API documentation (generated from docstrings)
- **[Basic Usage Guide](docs/usage/basic_usage_guide.md)** - Comprehensive step-by-step guide to using sprime with your data, including advanced usage patterns and testing

### Core Concepts

- **[Background and Concepts](docs/background_and_concepts.md)** - Introduction to qHTS, assays, S' metric, and key terminology
- **[Understanding 4PL Dose-Response Curves](docs/README_4PL_Dose_Response.md)** - Detailed explanation of the Hill equation model

### Technical Documentation

- **[Hill Curve Fitting Configuration](docs/usage/hill_curve_fitting_configuration.md)** - Technical details on curve fitting parameters and configuration options
- **[Development guide](docs/development.md)** - Dev setup, tests, pre-commit, API docs, CI (for contributors)

### Examples

- **[Demo Script](docs/usage/demo.py)** - Consolidated demo: S' only (raw -> Hill -> S'), pre-calc only (bypass raw data), and Delta S' (raw -> Hill -> S' -> delta, export tables)

## Requirements

### Runtime Dependencies

- **numpy** >= 1.20.0, < 3.0 - Numerical computing (compatible with numpy 1.x and 2.x)
- **scipy** >= 1.7.0 - Scientific computing and curve fitting

### Development Dependencies

Development dependencies are optional and only needed if you plan to contribute to sprime or run the test suite:

- **pytest** >= 7.0.0 - Testing framework for running the comprehensive test suite
  - Used to execute unit tests, integration tests, and edge case tests
  - Provides test discovery, fixtures, and assertion utilities
  - Required for: `pytest tests/` commands

- **pytest-cov** >= 4.0.0 - Coverage reporting plugin for pytest
  - Generates code coverage reports to identify untested code
  - Used with: `pytest --cov=src/sprime --cov-report=html`
  - Helps ensure comprehensive test coverage during development

- **pdoc3** - Generates API docs from docstrings into `pdoc_html/`
- **pre-commit** - Optional; rebuilds `pdoc_html` before each commit when `.py` files change

**Install development dependencies:**
```bash
pip install -e ".[dev]"
```

These dependencies are not required for normal usage of sprime - only for development and testing.

## Testing

sprime includes a comprehensive test suite. To run tests:

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run all tests
pytest tests/

# Run with coverage
pytest tests/ --cov=src/sprime --cov-report=html

# Run specific test files
pytest tests/test_sprime.py
pytest tests/test_hill_fitting.py
pytest tests/test_integration.py
```

See the [Basic Usage Guide](docs/usage/basic_usage_guide.md#running-the-test-suite) for more testing information.

### API docs (pdoc_html)

API docs live in `pdoc_html/` in the repo and are built from docstrings. The [API Reference](https://mocomakers.github.io/sprime/) on GitHub Pages is deployed from CI. Build locally: `python -m pdoc --html --output-dir pdoc_html --force sprime`. Optional: `pre-commit install` rebuilds and stages `pdoc_html` on commit when `.py` files change. See [Development guide](docs/development.md) for details.

## Project Status

This project is currently in **active development** (Alpha status). Features and API may change.

### Current Version

The package uses semantic versioning. Version is derived from **Git tags** (via [hatch-vcs](https://github.com/ofek/hatch-vcs)); `src/sprime/_version.py` is generated at build/install time and is not committed. Check the installed version:

```python
import sprime
print(sprime.__version__)
```

**Tagging releases:** Create a Git tag (e.g. `v0.1.0`) and push it. The version used in built packages and on PyPI comes from that tag.

```bash
git tag v0.1.0
git push origin v0.1.0
```

## Contributing

Contributions are welcome! We appreciate your help in making sprime better.

### How to Contribute

1. **Fork the repository** on GitHub
2. **Create a feature branch** (`git checkout -b feature/amazing-feature`)
3. **Set up dev environment** — venv, `pip install -e ".[dev]"`, `pytest tests/`, optionally `pre-commit install`. See [Development guide](docs/development.md).
4. **Make your changes** and add tests if applicable
5. **Run the test suite** — `pytest tests/`
6. **Commit your changes** (`git commit -m 'Add amazing feature'`). If you use pre-commit, the pdoc hook will rebuild and stage `pdoc_html` when you change `.py` files.
7. **Push to the branch** (`git push origin feature/amazing-feature`)
8. **Open a Pull Request**

### Development Guidelines

- Follow PEP 8 style guidelines
- Add tests for new features
- Update documentation as needed
- Ensure all tests pass before submitting

### Reporting Issues

If you find a bug or have a feature request, please open an issue on GitHub:
- [Issue Tracker](https://github.com/MoCoMakers/sprime/issues)

### Questions

For questions or support, please reach out via:
- [MoCo Makers Contact](https://www.mocomakers.com/contact/)
- Open a GitHub Discussion

## Repository

- **GitHub**: [https://github.com/MoCoMakers/sprime](https://github.com/MoCoMakers/sprime)
- **Issues**: [https://github.com/MoCoMakers/sprime/issues](https://github.com/MoCoMakers/sprime/issues)

## Citation

If you use sprime in your research, please cite:

### Generation 2 (S') Methodology

New citations introducing Generation 2 (S') are forthcoming. Please check back for updated citation information.

### Original S Metric

This library implements **Generation 2** of the S' methodology, which evolved from the original S metric described in:

> Zamora PO, Altay G, Santamaria U, Dwarshuis N, Donthi H, Moon CI, Bakalar D, Zamora M. Drug Responses in Plexiform Neurofibroma Type I (PNF1) Cell Lines Using High-Throughput Data and Combined Effectiveness and Potency. *Cancers (Basel)*. 2023 Dec 12;15(24):5811. doi: [10.3390/cancers15245811](https://doi.org/10.3390/cancers15245811). PMID: 38136356; PMCID: PMC10742026.

### Hill Curve Fitting Implementation

The Hill curve fitting implementation in sprime is inspired by the four parameter logistic regression work by Giuseppe Cardillo:

> Giuseppe Cardillo (2025). Four parameters logistic regression - There and back again (https://github.com/dnafinder/logistic4), GitHub. Retrieved December 24, 2025

## Acknowledgments

This library was developed by the **Computational Biology Working Group** at [DMV Petri Dish](https://compbio.dmvpetridish.com/) with support from the nonprofit [DMV Petri Dish](https://www.dmvpetridish.com/). 

We're also part of the DIY tech/science community at [MoCo Makers](https://www.mocomakers.com/).

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

Copyright (C) 2026 MoCo Maker Labs LLC

## Support

- **Documentation**: See the [docs](docs/) directory for comprehensive guides
- **Issues**: Report bugs or request features on [GitHub Issues](https://github.com/MoCoMakers/sprime/issues)
- **Contact**: Reach out via [MoCo Makers Contact](https://www.mocomakers.com/contact/)

## Related Projects

- [DMV Petri Dish Computational Biology](https://compbio.dmvpetridish.com/) - Computational biology research group
- [MoCo Makers](https://www.mocomakers.com/) - DIY tech/science community
