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

2. When you `git commit`, pre-commit runs the configured hooks. The **pdoc** hook:
   - Runs only if staged files include `.py` changes (`files: \.py$` in config).
   - Executes `python scripts/build_docs_precommit.py`.
   - That script runs `pdoc --html --output-dir pdoc_html --force sprime` then `git add pdoc_html`.
   - Your commit includes both your changes and updated `pdoc_html/`. Nothing touches PyPI.

**Config:** `.pre-commit-config.yaml` defines when (on commit, when `.py` changes) and what (`scripts/build_docs_precommit.py`).

**Script:** `scripts/build_docs_precommit.py` does the work: build docs, stage them. Used only by the pre-commit hook.

## Versioning

Version comes from **Git tags** via [hatch-vcs](https://github.com/ofek/hatch-vcs). `src/sprime/_version.py` is generated at build/install time and is **not** committed. To release:

```bash
git tag v0.1.0
git push origin v0.1.0
```

## CI

- **Deploy API docs:** `.github/workflows/deploy-api-docs.yml` runs on push to `main` or tags `v*`. It builds pdoc into `build/sprime`, uploads to GitHub Pages. The live API Reference is always the latest deploy. No versioned doc paths.
