"""
Golden checks for committed CSVs under docs/usage/.

Expected S' values are **library output** after load -> process. The
``demo_precontrol_normalized_precalc.csv`` ``S'`` column matches sprime's
``asinh((Zero - Inf) / AC50)`` recomputation (rounded in the CSV; process may differ
in the last digits).

Raw vehicle demos use percent nucleus + DMSO 35.3 from ``tests/fixtures/SPrime_variation_reference.csv``.
Default ``load`` + ``process`` applies ``pipeline_asymptote_normalized`` (``normalize_to_max_value`` + x100) before fitting.
S' = ``asinh((Zero-Inf)/EC50)`` (signed).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from sprime import SPrime

_USAGE_DIR = Path(__file__).resolve().parent.parent / "docs" / "usage"

# Locked against sprime + scipy in CI (2025-03); rtol catches unintended pipeline changes.
_RTOL = 1e-4


def _s_prime_by_cell_line(path: Path, **load_kw):
    load_kw.setdefault("response_normalization", "asymptote_normalized")
    raw, _ = SPrime.load(path, **load_kw)
    screening, _ = SPrime.process(raw, allow_overwrite_precalc_params=True)
    return {p.cell_line.name: p.s_prime for p in screening.profiles}


class TestPrecontrolNormalizedDemos:
    """``skip_control_response_normalization=True`` - empty Control_Response, ipNF-style scale."""

    @pytest.fixture
    def skip(self):
        return {
            "skip_control_response_normalization": True,
            "response_normalization": "asymptote_normalized",
        }

    def test_demo_precontrol_normalized_s_prime_csv(self, skip):
        sp_by_cl = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_precontrol_normalized_s_prime.csv",
            values_as="columns",
            **skip,
        )
        assert sp_by_cl["ipNF96.11C"] == pytest.approx(-18.848790037581885, rel=_RTOL)

    def test_demo_precontrol_normalized_raw_list_csv(self, skip):
        sp_by_cl = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_precontrol_normalized_raw_list.csv",
            values_as="list",
            **skip,
        )
        assert sp_by_cl["ipNF96.11C"] == pytest.approx(-18.84880324391443, rel=_RTOL)

    def test_demo_precontrol_normalized_raw_list_matches_columns(self, skip):
        """List vs column layout should yield the same S' for the same biology."""
        cols = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_precontrol_normalized_s_prime.csv",
            values_as="columns",
            **skip,
        )
        lst = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_precontrol_normalized_raw_list.csv",
            values_as="list",
            **skip,
        )
        assert cols["ipNF96.11C"] == pytest.approx(lst["ipNF96.11C"], abs=1e-3)

    def test_demo_precontrol_normalized_delta_csv(self, skip):
        sp_by_cl = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_precontrol_normalized_delta.csv",
            values_as="columns",
            **skip,
        )
        assert sp_by_cl["ipNF96.11C"] == pytest.approx(-18.848790037581885, rel=_RTOL)
        assert sp_by_cl["ipnf05.5 mc"] == pytest.approx(-18.858763650417078, rel=_RTOL)

        raw, _ = SPrime.load(
            _USAGE_DIR / "demo_precontrol_normalized_delta.csv",
            values_as="columns",
            **skip,
        )
        screening, _ = SPrime.process(raw, allow_overwrite_precalc_params=True)
        delta = screening.calculate_delta_s_prime(
            reference_cell_lines="ipnf05.5 mc",
            test_cell_lines=["ipNF96.11C"],
        )
        row = delta["ipnf05.5 mc"][0]
        assert row["delta_s_prime"] == pytest.approx(-0.009973612835192824, rel=_RTOL)

    def test_demo_precontrol_normalized_precalc_csv(self):
        """Path B only; S' matches CSV reference (sprime recomputation from Hill params)."""
        sp_by_cl = _s_prime_by_cell_line(_USAGE_DIR / "demo_precontrol_normalized_precalc.csv")
        assert sp_by_cl["ipNF96.11C"] == pytest.approx(-5.044950504499101, rel=_RTOL)
        assert sp_by_cl["ipnf05.5 mc"] == pytest.approx(-5.046821041, rel=_RTOL)
        assert sp_by_cl["ipnf02.3"] == pytest.approx(-5.032656731, rel=_RTOL)


class TestRawVehicleControlDemos:
    """Default strict load + asymptote-normalized pipeline; S' = asinh((Zero-Inf)/EC50)."""

    def test_demo_raw_vehicle_control_s_prime_csv(self):
        sp_by_cl = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_raw_vehicle_control_s_prime.csv",
            values_as="columns",
        )
        assert sp_by_cl["NF2-/-"] == pytest.approx(6.867631094104203, rel=_RTOL)

    def test_demo_raw_vehicle_control_raw_list_csv(self):
        sp_by_cl = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_raw_vehicle_control_raw_list.csv",
            values_as="list",
        )
        assert sp_by_cl["NF2-/-"] == pytest.approx(6.867631094104203, rel=_RTOL)

    def test_demo_raw_vehicle_control_raw_list_matches_columns(self):
        cols = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_raw_vehicle_control_s_prime.csv",
            values_as="columns",
        )
        lst = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_raw_vehicle_control_raw_list.csv",
            values_as="list",
        )
        assert cols["NF2-/-"] == pytest.approx(lst["NF2-/-"], abs=1e-6)

    def test_demo_raw_vehicle_control_delta_csv(self):
        sp_by_cl = _s_prime_by_cell_line(
            _USAGE_DIR / "demo_raw_vehicle_control_delta.csv",
            values_as="columns",
        )
        assert sp_by_cl["NF2-/- Ref"] == pytest.approx(-3.1764092497421257, rel=_RTOL)
        assert sp_by_cl["NF2-/- Test_1"] == pytest.approx(-6.143367466705312, rel=_RTOL)
        assert sp_by_cl["NF2-/- Test_2"] == pytest.approx(-4.720064185865709, rel=_RTOL)

        raw, _ = SPrime.load(
            _USAGE_DIR / "demo_raw_vehicle_control_delta.csv",
            values_as="columns",
            response_normalization="asymptote_normalized",
        )
        screening, _ = SPrime.process(raw)
        delta = screening.calculate_delta_s_prime(
            reference_cell_lines="NF2-/- Ref",
            test_cell_lines=["NF2-/- Test_1", "NF2-/- Test_2"],
        )
        rows = delta["NF2-/- Ref"]
        by_test = {r["test_cell_line"]: r["delta_s_prime"] for r in rows}
        assert by_test["NF2-/- Test_1"] == pytest.approx(2.9669582169631865, rel=_RTOL)
        assert by_test["NF2-/- Test_2"] == pytest.approx(1.5436549361235836, rel=_RTOL)
