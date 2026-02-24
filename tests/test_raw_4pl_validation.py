"""
Tests validating sprime's 4PL fitting and S' calculation against Paul's
reference Excel workbook.

Reference data source:
    data/working 1___ctf-ucf-cenix-completecompounddata-all passes in one
    tab_ to release_final_cleaned (3).xlsx
    Sheet: 'Sprime_delta_Sprime'

Three representative dose-response profiles are extracted:
    1. 17-AAG in NF2-/- at 48h  (rows 7-14)
    2. AR42  in NF2-/- at 48h  (rows 19-26)
    3. AR42  in wt     at 48h  (rows 31-37)

Each profile is provided with both:
    (a) raw responses and DMSO control values from the Excel sheet, and
    (b) the already-normalised responses (scaled 0-100) that Paul
        computed before 4PL fitting.

The expected 4PL parameters (A, B, C, D) and S' values come from the
same Excel sheet.
"""

import math
import pytest

from sprime.hill_fitting import fit_hill_curve
from sprime.sprime import (
    HillCurveParams,
    DoseResponseProfile,
    Compound,
    CellLine,
    Assay,
    calculate_s_prime_from_params,
    normalize_responses,
)


# ============================================================================
# Reference data extracted from Paul's Excel
# ============================================================================

# 17-AAG in NF2-/- at 48 hours (rows 7-14, sheet Sprime_delta_Sprime)
AAG_17_NF2 = {
    "compound": "17-AAG",
    "cell_line": "NF2-/-",
    "concentrations_uM": [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0],
    # Raw responses (col K: TripMean_number_of_nuclei equivalent)
    "raw_responses": [
        41.2416514672222,
        40.8554400777777,
        38.7785148277777,
        38.3814239355555,
        6.63547943372222,
        0.936100936333333,
        0.0,
        0.0,
    ],
    # DMSO control (row 15, col K): single replicate
    "dmso_control": 35.3,
    # Already-normalised responses (col O: x 100 scalar)
    "responses_normalized": [
        100.0,
        99.06354043617421,
        94.02755090590361,
        93.0647114508983,
        16.089267033828477,
        2.269794984028033,
        0.0,
        0.0,
    ],
    # Paul's 4PL: y = D + (A-D)/(1+(x/C)^B)
    "paul_A": 97.736,       # upper asymptote (response at low conc)
    "paul_B": 4.2145,       # Hill coefficient (positive convention)
    "paul_C_ec50": 0.202,   # EC50 (µM)
    "paul_D": 0.7084,       # lower asymptote (response at high conc)
    "paul_emax": 97.0276,   # A - D
    "paul_s_prime": 6.867631319810275,
}

# AR42 in NF2-/- at 48 hours (rows 19-26)
AR42_NF2 = {
    "compound": "AR42",
    "cell_line": "NF2-/-",
    "concentrations_uM": [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0],
    "raw_responses": [
        32.3333333333333,
        30.9444444444444,
        26.1111111111111,
        17.9444444444444,
        8.5,
        0.0,
        0.0,
        0.0,
    ],
    "dmso_control": 35.3,
    "responses_normalized": [
        100.0,
        95.70446735395186,
        80.75601374570451,
        55.498281786941504,
        26.28865979381446,
        0.0,
        0.0,
        0.0,
    ],
    "paul_A": 100.677,
    "paul_B": 1.1465,
    "paul_C_ec50": 0.1216,
    "paul_D": -3.2032,
    "paul_emax": 103.8802,
    "paul_s_prime": 7.443404144664863,
}

# AR42 in wt (wild type) at 48 hours (rows 31-37)
AR42_WT = {
    "compound": "AR42",
    "cell_line": "wt",
    "concentrations_uM": [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0],
    "raw_responses": [
        2.44444444444444,
        2.44444444444444,
        1.94444444444444,
        0.555555555555555,
        0.0,
        0.0,
        0.0,
    ],
    # DMSO average of two replicates: (3.16666666666666 + 4.375) / 2
    "dmso_control": 3.77083333333333,
    "responses_normalized": [
        100.0,
        100.0,
        79.5454545454545,
        22.727272727272744,
        0.0,
        0.0,
        0.0,
    ],
    "paul_A": 101.0881,
    "paul_B": 2.1603,
    "paul_C_ec50": 0.056,
    "paul_D": -0.5898,
    "paul_emax": 101.6779,
    "paul_s_prime": 8.197360818279554,
}

DATASETS = [AAG_17_NF2, AR42_NF2, AR42_WT]


def _dataset_id(d):
    """Readable test-ID for parametrize."""
    return f"{d['compound']}_{d['cell_line']}"


# ============================================================================
# Normalization tests
# ============================================================================


class TestDMSONormalization:
    """Verify that normalize_responses reproduces Paul's Excel preprocessing."""

    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_normalize_matches_excel(self, dataset):
        """normalize_responses(raw, dmso) should match Paul's scaled column."""
        normalized = normalize_responses(
            dataset["raw_responses"],
            dataset["dmso_control"],
        )
        for i, (got, expected) in enumerate(
            zip(normalized, dataset["responses_normalized"])
        ):
            assert got == pytest.approx(expected, abs=1e-8), (
                f"Mismatch at index {i}: got {got}, expected {expected}"
            )

    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_normalize_length_preserved(self, dataset):
        """Output length must match input length."""
        normalized = normalize_responses(
            dataset["raw_responses"],
            dataset["dmso_control"],
        )
        assert len(normalized) == len(dataset["raw_responses"])

    def test_normalize_zero_control_skips_division(self):
        """control_value=0 skips the control division step (no error)."""
        result = normalize_responses([10, 5, 0], 0.0)
        # Skips control division; max-normalises raw values, then scales
        assert result[0] == pytest.approx(100.0)
        assert result[1] == pytest.approx(50.0)
        assert result[2] == pytest.approx(0.0)

    def test_normalize_negative_control_returns_zeros(self):
        """Negative control_value flips signs; max ratio is 0 → all zeros."""
        result = normalize_responses([10, 5, 0], -1.0)
        # Ratios become [-10, -5, 0]; max(ratios) = 0 → returns zeros
        assert result == [0.0, 0.0, 0.0]

    def test_normalize_all_zero_responses(self):
        """All-zero responses return all zeros (no error)."""
        result = normalize_responses([0, 0, 0], 1.0)
        assert result == [0.0, 0.0, 0.0]

    def test_normalize_custom_scale(self):
        """Scaling factor of 1.0 gives [0, 1] range."""
        normalized = normalize_responses([10.0, 5.0, 0.0], 10.0, scale=1.0)
        assert normalized[0] == pytest.approx(1.0)
        assert normalized[1] == pytest.approx(0.5)
        assert normalized[2] == pytest.approx(0.0)

    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_normalize_then_assign_to_profile(self, dataset):
        """normalize_responses → assign to profile.responses should work."""
        normalized = normalize_responses(
            dataset["raw_responses"],
            dataset["dmso_control"],
        )
        profile = DoseResponseProfile(
            compound=Compound(name=dataset["compound"], drug_id=dataset["compound"]),
            cell_line=CellLine(name=dataset["cell_line"]),
            assay=Assay(name="test"),
            concentrations=dataset["concentrations_uM"],
            responses=normalized,
            concentration_units="microM",
        )
        for got, expected in zip(profile.responses, dataset["responses_normalized"]):
            assert got == pytest.approx(expected, abs=1e-8)


# ============================================================================
# 4PL fitting + S' validation (using pre-normalised data)
# ============================================================================


class TestPaulReferenceRawTo4PL:
    """Validate sprime's 4PL fitting against Paul's Excel reference values."""

    # ------------------------------------------------------------------
    # EC50
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_fitted_ec50_matches(self, dataset):
        """Fitted EC50 should match Paul's value within 15 % relative."""
        params = fit_hill_curve(
            dataset["concentrations_uM"],
            dataset["responses_normalized"],
        )
        assert params.ec50 == pytest.approx(dataset["paul_C_ec50"], rel=0.15)

    # ------------------------------------------------------------------
    # Emax (|Upper - Lower|)
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_fitted_emax_matches(self, dataset):
        """Absolute amplitude |Upper - Lower| should match Paul's Emax within 5 %."""
        params = fit_hill_curve(
            dataset["concentrations_uM"],
            dataset["responses_normalized"],
        )
        fitted_emax = abs(params.upper - params.lower)
        assert fitted_emax == pytest.approx(dataset["paul_emax"], rel=0.05)

    # ------------------------------------------------------------------
    # Hill coefficient (absolute value; sign depends on parameterisation)
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_fitted_hill_coefficient_magnitude_matches(self, dataset):
        """|Hill coefficient| should match Paul's B within 15 %."""
        params = fit_hill_curve(
            dataset["concentrations_uM"],
            dataset["responses_normalized"],
        )
        assert abs(params.hill_coefficient) == pytest.approx(
            dataset["paul_B"], rel=0.15
        )

    # ------------------------------------------------------------------
    # R² (goodness of fit)
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_fit_quality(self, dataset):
        """R² should be > 0.99 for these well-behaved profiles."""
        params = fit_hill_curve(
            dataset["concentrations_uM"],
            dataset["responses_normalized"],
        )
        assert params.r_squared is not None
        assert params.r_squared > 0.99

    # ------------------------------------------------------------------
    # S' (the primary output metric)
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_s_prime_matches(self, dataset):
        """S' = asinh(|Upper-Lower| / EC50) should match Paul's value within 1 %."""
        params = fit_hill_curve(
            dataset["concentrations_uM"],
            dataset["responses_normalized"],
        )
        # Use abs(amplitude) so sign is always positive, matching Paul
        s_prime = math.asinh(abs(params.upper - params.lower) / params.ec50)
        assert s_prime > 0, "S' should always be positive"
        assert s_prime == pytest.approx(dataset["paul_s_prime"], rel=0.01)

    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_s_prime_positive_via_profile(self, dataset):
        """DoseResponseProfile.calculate_s_prime() should return a positive S'."""
        profile = DoseResponseProfile(
            compound=Compound(name=dataset["compound"], drug_id=dataset["compound"]),
            cell_line=CellLine(name=dataset["cell_line"]),
            assay=Assay(name="test"),
            concentrations=dataset["concentrations_uM"],
            responses=dataset["responses_normalized"],
            concentration_units="microM",
        )
        s_prime = profile.fit_and_calculate_s_prime()
        assert s_prime > 0, "S' via DoseResponseProfile should always be positive"
        assert s_prime == pytest.approx(dataset["paul_s_prime"], rel=0.01)

    # ------------------------------------------------------------------
    # Delta S' end-to-end (AR42: wt vs NF2-/-)
    # ------------------------------------------------------------------
    def test_delta_s_prime_ar42(self):
        """Delta S' = S'(wt) - S'(NF2-/-) should match Paul's value."""
        expected_delta = 0.7539566736146908  # from Excel row 44

        profiles = []
        for ds in [AR42_NF2, AR42_WT]:
            profile = DoseResponseProfile(
                compound=Compound(name=ds["compound"], drug_id=ds["compound"]),
                cell_line=CellLine(name=ds["cell_line"]),
                assay=Assay(name="test"),
                concentrations=ds["concentrations_uM"],
                responses=ds["responses_normalized"],
                concentration_units="microM",
            )
            profile.fit_and_calculate_s_prime()
            profiles.append(profile)

        s_prime_nf2 = next(p for p in profiles if p.cell_line.name == "NF2-/-").s_prime
        s_prime_wt = next(p for p in profiles if p.cell_line.name == "wt").s_prime

        delta = s_prime_wt - s_prime_nf2
        assert delta == pytest.approx(expected_delta, abs=0.05)

    # ------------------------------------------------------------------
    # Standalone calculate_s_prime_from_params
    # ------------------------------------------------------------------
    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_calculate_s_prime_from_params_matches(self, dataset):
        """calculate_s_prime_from_params should match Paul's S' when fed his params."""
        s_prime = calculate_s_prime_from_params(
            ac50=dataset["paul_C_ec50"],
            upper=dataset["paul_A"],
            lower=dataset["paul_D"],
        )
        assert s_prime > 0
        assert s_prime == pytest.approx(dataset["paul_s_prime"], rel=0.001)

    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_calculate_s_prime_from_params_swapped(self, dataset):
        """calculate_s_prime_from_params with upper/lower swapped still gives positive S'."""
        s_prime = calculate_s_prime_from_params(
            ac50=dataset["paul_C_ec50"],
            upper=dataset["paul_D"],    # swapped
            lower=dataset["paul_A"],    # swapped
        )
        assert s_prime > 0
        assert s_prime == pytest.approx(dataset["paul_s_prime"], rel=0.001)


# ============================================================================
# End-to-end: raw → normalize → fit → S' → ΔS'
# ============================================================================


class TestRawToSPrimeEndToEnd:
    """
    Full pipeline test: feed raw (unnormalised) responses, apply
    normalize_responses, then fit and compute S' / ΔS'.

    This verifies that the library can reproduce Paul's Excel results
    starting from the raw instrument readouts.
    """

    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_raw_normalize_fit_s_prime(self, dataset):
        """Raw → normalise → fit → S' should match Paul's reference."""
        # 1. Normalise
        normalized = normalize_responses(
            dataset["raw_responses"],
            dataset["dmso_control"],
        )

        # 2. Build profile with normalised data
        profile = DoseResponseProfile(
            compound=Compound(name=dataset["compound"], drug_id=dataset["compound"]),
            cell_line=CellLine(name=dataset["cell_line"]),
            assay=Assay(name="test"),
            concentrations=dataset["concentrations_uM"],
            responses=normalized,
            concentration_units="microM",
        )

        # 3. Fit + S'
        s_prime = profile.fit_and_calculate_s_prime()

        # 4. Assertions
        assert s_prime > 0
        assert s_prime == pytest.approx(dataset["paul_s_prime"], rel=0.01)
        assert profile.hill_params.ec50 == pytest.approx(
            dataset["paul_C_ec50"], rel=0.15
        )

    @pytest.mark.parametrize("dataset", DATASETS, ids=_dataset_id)
    def test_normalize_assign_then_fit(self, dataset):
        """normalize_responses → profile.responses → fit_and_calculate_s_prime()."""
        normalized = normalize_responses(
            dataset["raw_responses"],
            dataset["dmso_control"],
        )
        profile = DoseResponseProfile(
            compound=Compound(name=dataset["compound"], drug_id=dataset["compound"]),
            cell_line=CellLine(name=dataset["cell_line"]),
            assay=Assay(name="test"),
            concentrations=dataset["concentrations_uM"],
            responses=normalized,
            concentration_units="microM",
        )

        # Fit + S'
        s_prime = profile.fit_and_calculate_s_prime()
        assert s_prime > 0
        assert s_prime == pytest.approx(dataset["paul_s_prime"], rel=0.01)

    def test_raw_to_delta_s_prime_ar42(self):
        """Full pipeline: raw → normalise → fit → ΔS'(wt – NF2-/-) for AR42."""
        expected_delta = 0.7539566736146908  # from Excel row 44

        profiles = []
        for ds in [AR42_NF2, AR42_WT]:
            normalized = normalize_responses(
                ds["raw_responses"],
                ds["dmso_control"],
            )
            profile = DoseResponseProfile(
                compound=Compound(name=ds["compound"], drug_id=ds["compound"]),
                cell_line=CellLine(name=ds["cell_line"]),
                assay=Assay(name="test"),
                concentrations=ds["concentrations_uM"],
                responses=normalized,
                concentration_units="microM",
            )
            profile.fit_and_calculate_s_prime()
            profiles.append(profile)

        s_prime_nf2 = next(p for p in profiles if p.cell_line.name == "NF2-/-").s_prime
        s_prime_wt = next(p for p in profiles if p.cell_line.name == "wt").s_prime

        delta = s_prime_wt - s_prime_nf2
        assert delta == pytest.approx(expected_delta, abs=0.05)
