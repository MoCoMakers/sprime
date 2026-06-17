#!/usr/bin/env bash
set -e
python -V
pip install -e ".[dev]"
pytest -q



python scripts/excel_probe.py
