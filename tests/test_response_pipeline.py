"""Tests for response preprocessing aligned with SPrime_variation_reference.csv."""

import csv
from pathlib import Path

import pytest

from sprime.response_pipeline import (
    S_PRIME_RESPONSE_SCALE_FACTOR,
    normalize_to_max_value,
    pipeline_asymptote_normalized,
    pipeline_response_scale,
    ratios_to_control,
    scale_responses,
)

# Values from tests/fixtures/SPrime_variation_reference.csv (update literals/tolerances if the fixture changes).
_REFERENCE_RAW_17_AAG = [41.24, 40.86, 38.78, 38.38, 6.64, 0.94, 0.00, 0.00]
_REFERENCE_DMSO = 35.3
_FIXTURE_CSV = Path(__file__).resolve().parent / "fixtures" / "SPrime_variation_reference.csv"


def test_reference_csv_row2_intermediates_17_aag():
    """First data row of SPrime_variation_reference.csv (17-AAG, 0.003 microM)."""
    raw = _REFERENCE_RAW_17_AAG[0]
    dmso = _REFERENCE_DMSO
    ratio = ratios_to_control([raw], dmso)[0]
    assert abs(ratio - raw / dmso) < 1e-12
    # Tabulated ratio uses more digits than percent/control imply; stay close to the sheet.
    assert abs(ratio - 1.168318738) < 5e-4

    ratios = ratios_to_control(_REFERENCE_RAW_17_AAG, dmso)
    unity = normalize_to_max_value(ratios)
    assert abs(unity[0] - 1.0) < 1e-9
    assert abs(unity[1] - 0.990635404) < 2e-4

    scaled = scale_responses(unity)
    assert abs(scaled[0] - 100.0) < 1e-6
    assert abs(scaled[1] - 99.06) < 0.05

    assert S_PRIME_RESPONSE_SCALE_FACTOR == 100.0


def test_pipeline_response_scale_matches_direct_x100():
    """Response-scale path equals 100 * (raw / control)."""
    out = pipeline_response_scale(_REFERENCE_RAW_17_AAG, _REFERENCE_DMSO)
    for r, o in zip(_REFERENCE_RAW_17_AAG, out):
        if r == 0.0:
            assert abs(o) < 1e-12
        else:
            assert abs(o - 100.0 * r / _REFERENCE_DMSO) < 1e-9


def test_ratios_to_control_zero_raises():
    with pytest.raises(ValueError, match="non-zero"):
        ratios_to_control([1.0, 2.0], 0.0)


def test_normalize_to_max_value_zero_max_raises():
    with pytest.raises(ValueError, match="maximum value is zero"):
        normalize_to_max_value([0.0, 0.0])


def test_normalize_to_max_value_empty_raises():
    with pytest.raises(ValueError, match="non-empty"):
        normalize_to_max_value([])


def test_pipeline_asymptote_normalized_matches_golden_path():
    out = pipeline_asymptote_normalized(_REFERENCE_RAW_17_AAG, _REFERENCE_DMSO)
    assert abs(out[0] - 100.0) < 1e-6
    assert abs(out[1] - 99.06) < 0.02


def test_reference_csv_internal_consistency():
    """
    Parse the reference CSV (human-maintained): Nonnormalized x100 ~= 100 x second ratio column.

    Percent nucleus and ratio columns may disagree slightly from float(raw/control) because
    the sheet rounds display values; we do not assert bit-exact pipeline equality to every column.
    """
    path = _FIXTURE_CSV
    with path.open(encoding="latin-1", newline="") as f:
        rows = list(csv.reader(f))
    header, body = rows[0], rows[1:]

    def col_index(name: str) -> int:
        try:
            return header.index(name)
        except ValueError as e:
            raise AssertionError(f"column {name!r} not in header") from e

    i_r2 = len(header) - 1 - header[::-1].index("test/control ratio")
    i_nn = col_index("Nonnormalized, x100 scale")

    for row in body:
        if not row or not row[0].strip() or row[0].strip().upper() == "DMSO":
            continue
        r2_s = row[i_r2].strip()
        nn_s = row[i_nn].strip()
        if r2_s == "" or nn_s == "":
            continue
        ratio_nn = float(r2_s)
        nn = float(nn_s)
        assert abs(nn - 100.0 * ratio_nn) < 1e-5, (row[0], nn, ratio_nn)
