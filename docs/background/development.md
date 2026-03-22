# Developer setup

Details new contributors need to run tests, build API docs, and use pre-commit.

## Environment

```bash
git clone https://github.com/MoCoMakers/sprime.git
cd sprime
python -m venv venv
# Windows: venv\Scripts\activate
# macOS/Linux: source venv/bin/activate
pip install -e ".[dev]"
```

## Product / pipeline documentation

Behavioral branches for S' (pre-calculated vs raw, **Control_Response** / DMSO-relative ratios, asymptote-normalized vs response-scale, x100 convention, and **`skip_control_response_normalization`** on loaders) are documented in **[`s_prime_derivation_pipeline.md`](s_prime_derivation_pipeline.md)** (under **`docs/background/`**). The x100-aligned helpers live in **`sprime.response_pipeline`**. Pre-calculated overwrite behavior uses **`allow_overwrite_precalc_params`** (replaces the old `allow_overwrite_hill_coefficients` name). Update docs when API or validation rules change.

**Variation reference:** `tests/fixtures/SPrime_variation_reference.csv` is the committed copy used by `tests/test_response_pipeline.py`. Update it when reference numbers change. **Private scratch** under `docs/` is listed in **`.gitignore`** and is not published--do not link to those paths from **README** or user-facing **docs/**.

## Linting and formatting (Ruff)

This project uses **[Ruff](https://docs.astral.sh/ruff/)** for style and static checks. Rules and line length are defined in **`pyproject.toml`** (`[tool.ruff]` / `[tool.ruff.lint]` / `[tool.ruff.format]`).

After `pip install -e ".[dev]"`, run Ruff from the **repository root** so lint/format match **pre-commit** (which can include **`docs/`** notebooks, not only `src/` and `tests/`):

```bash
ruff check .                  # recommended: whole tree (respects .gitignore, e.g. skips venv/)
ruff check . --fix            # apply safe auto-fixes
ruff format .                 # formatter (optional; not run by pre-commit by default)
```

Narrower scope when debugging: `ruff check src tests`.

## Tests

```bash
pytest tests/
pytest tests/ --cov=src/sprime --cov-report=html   # with coverage
```

## API docs (`pdoc_html`)

API docs live in `pdoc_html/` in the repo and are built from docstrings. The live [API Reference](https://mocomakers.github.io/sprime/) is deployed from CI (GitHub Pages). You need to do `pip install pdoc3`

**Build locally:**

```bash
python -m pdoc --html --output-dir pdoc_html --force sprime
```

**Update on commit:** A pre-commit hook can rebuild and stage `pdoc_html` automatically when you change `.py` files. See [Pre-commit](#pre-commit) below.

## Pre-commit

Git does **not** read `.pre-commit-config.yaml`. The [pre-commit](https://pre-commit.com/) tool does. You install a Git hook that runs pre-commit on each `git commit`.

1. **Install the hook** (once per clone):

   ```bash
   pre-commit install
   ```

2. When you `git commit`, pre-commit runs the configured hooks **before** Git finishes the commit:

   - **Ruff** -- runs `ruff check` on **staged** paths (see `.pre-commit-config.yaml`). That can include e.g. **`docs/**/*.ipynb`**, not only `src/` and `tests/`—so **`ruff check .`** from the repo root before committing aligns with what the hook may see. If Ruff reports any problem, the hook **exits with an error** and **Git aborts the commit** (your staged changes stay staged; fix or `ruff check . --fix`, then try again). Invalid Python is reported like other check failures.
   - **pdoc** -- runs only when staged files include `.py` changes (`files: \.py$`). It runs `python scripts/build_docs_precommit.py`, which builds `pdoc_html/` and stages it so the commit can include updated API HTML.

   To **skip all hooks** (not recommended): `git commit --no-verify`.

**Config:** `.pre-commit-config.yaml` lists hooks and versions.

**Script:** `scripts/build_docs_precommit.py` is used only by the pdoc hook.

## Versioning

Version comes from **Git tags** via [hatch-vcs](https://github.com/ofek/hatch-vcs). `src/sprime/_version.py` is generated at build/install time and is **not** committed. To release:

```bash
git tag v0.1.0
git push origin v0.1.0
```

## CI

- **Deploy API docs:** `.github/workflows/deploy-api-docs.yml` runs on push to `main` or tags `v*`. It builds pdoc into `build/sprime`, uploads to GitHub Pages. The live API Reference is always the latest deploy. No versioned doc paths.
