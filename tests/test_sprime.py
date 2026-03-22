"""
Tests for sprime core module.

Tests cover:
- RawDataset loading and processing
- ScreeningDataset operations
- SPrime API methods
- DoseResponseProfile operations
- Utility functions
- Edge cases and error handling
"""

import csv
import math
import tempfile
from pathlib import Path
from typing import Dict, List

import pytest

from sprime import (
    Assay,
    CellLine,
    Compound,
    DoseResponseProfile,
    HillCurveParams,
    RawDataset,
    ScreeningDataset,
    SPrime,
    calculate_delta_s_prime,
    calculate_s_prime_from_params,
    convert_to_micromolar,
    fit_hill_from_raw_data,
    get_s_prime_from_data,
    get_s_primes_from_file,
)
from sprime.response_pipeline import (
    pipeline_asymptote_normalized,
    pipeline_response_scale,
)


class TestCompound:
    """Tests for Compound value object."""

    def test_compound_creation(self):
        compound = Compound(name="Test Drug", drug_id="TEST001")
        assert compound.name == "Test Drug"
        assert compound.drug_id == "TEST001"
        assert compound.pubchem_sid is None
        assert compound.smiles is None

    def test_compound_with_optional_fields(self):
        compound = Compound(name="Test Drug", drug_id="TEST001", pubchem_sid="12345", smiles="CCO")
        assert compound.pubchem_sid == "12345"
        assert compound.smiles == "CCO"


class TestCellLine:
    """Tests for CellLine value object."""

    def test_cell_line_creation(self):
        cell_line = CellLine(name="Test Cell Line")
        assert cell_line.name == "Test Cell Line"
        assert cell_line.ref_id is None

    def test_cell_line_with_ref_id(self):
        cell_line = CellLine(name="Test Cell Line", ref_id="REF001")
        assert cell_line.ref_id == "REF001"


class TestAssay:
    """Tests for Assay value object."""

    def test_assay_creation(self):
        assay = Assay(name="Test Assay")
        assert assay.name == "Test Assay"
        assert assay.description is None
        assert assay.screen_id is None


class TestHillCurveParams:
    """Tests for HillCurveParams."""

    def test_hill_curve_params_creation(self):
        params = HillCurveParams(
            ec50=10.0,
            zero_asymptote=0.0,
            inf_asymptote=100.0,
            steepness_coefficient=1.5,
            r_squared=0.95,
        )
        assert params.ec50 == 10.0
        assert params.zero_asymptote == 0.0
        assert params.inf_asymptote == 100.0
        # Hill numerator in hill_equation uses (zero_asymptote - inf_asymptote); S' uses same order.
        assert params.inf_asymptote - params.zero_asymptote == 100.0
        assert params.zero_asymptote - params.inf_asymptote == -100.0

    def test_hill_curve_params_optional_fields(self):
        params = HillCurveParams(ec50=10.0, zero_asymptote=0.0, inf_asymptote=100.0)
        assert params.steepness_coefficient is None
        assert params.r_squared is None


class TestDoseResponseProfile:
    """Tests for DoseResponseProfile."""

    @pytest.fixture
    def profile(self):
        compound = Compound(name="Test Drug", drug_id="TEST001")
        cell_line = CellLine(name="Test Cell Line")
        assay = Assay(name="Test Assay")
        return DoseResponseProfile(
            compound=compound,
            cell_line=cell_line,
            assay=assay,
            concentrations=[0.1, 1, 10, 100, 1000],
            responses=[5, 10, 50, 90, 95],
        )

    def test_profile_creation(self, profile):
        assert profile.compound.name == "Test Drug"
        assert profile.cell_line.name == "Test Cell Line"
        assert len(profile.concentrations) == 5
        assert profile.s_prime is None
        assert profile.hill_params is None

    def test_fit_hill_curve(self, profile):
        params = profile.fit_hill_curve()
        assert isinstance(params, HillCurveParams)
        assert params.ec50 > 0
        assert params.r_squared is not None
        assert 0 <= params.r_squared <= 1

    def test_calculate_s_prime(self, profile):
        profile.fit_hill_curve()
        s_prime = profile.calculate_s_prime()
        assert s_prime is not None
        assert isinstance(s_prime, float)
        assert math.isfinite(s_prime)  # S' may be positive, negative, or ~zero depending on fit

    def test_fit_and_calculate_s_prime(self, profile):
        s_prime = profile.fit_and_calculate_s_prime()
        assert s_prime is not None
        assert profile.hill_params is not None
        assert profile.s_prime == s_prime

    def test_fit_hill_curve_without_data(self):
        compound = Compound(name="Test", drug_id="TEST")
        cell_line = CellLine(name="Test")
        assay = Assay(name="Test")
        profile = DoseResponseProfile(compound=compound, cell_line=cell_line, assay=assay)
        with pytest.raises(ValueError, match="Need raw data"):
            profile.fit_hill_curve()

    def test_calculate_s_prime_without_fit(self, profile):
        with pytest.raises(ValueError, match="Must fit Hill curve"):
            profile.calculate_s_prime()


class TestRawDataset:
    """Tests for RawDataset."""

    def create_test_csv(self, rows: List[Dict]) -> Path:
        """Helper to create temporary CSV file."""
        temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".csv", newline="")
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())

        writer = csv.DictWriter(temp_file, fieldnames=sorted(fieldnames))
        writer.writeheader()
        writer.writerows(rows)
        temp_file.close()
        return Path(temp_file.name)

    def test_load_from_file_with_raw_data(self):
        rows = [
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 1",
                "Concentration_Units": "microM",
                "Data0": "10",
                "Data1": "20",
                "Data2": "50",
                "Data3": "90",
                "Conc0": "0.1",
                "Conc1": "1",
                "Conc2": "10",
                "Conc3": "100",
            }
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, report = RawDataset.load_from_file(
                csv_file,
                skip_control_response_normalization=True,
                response_normalization="asymptote_normalized",
            )
            assert len(raw_data) == 1
            profile = list(raw_data.profiles)[0]
            assert profile.compound.name == "Drug A"
            assert len(profile.concentrations) == 4
        finally:
            csv_file.unlink()

    def test_load_from_file_with_precalculated_params(self):
        rows = [
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 1",
                "AC50": "10.0",
                "Upper": "100.0",
                "Lower": "0.0",
                "Hill_Slope": "1.5",
                "r2": "0.95",
            }
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, report = RawDataset.load_from_file(
                csv_file,
                skip_control_response_normalization=True,
                response_normalization="asymptote_normalized",
            )
            assert len(raw_data) == 1
            profile = list(raw_data.profiles)[0]
            assert profile.hill_params is not None
            assert profile.hill_params.ec50 == 10.0
        finally:
            csv_file.unlink()

    def test_rows_literal_empty_is_null(self):
        """Rows are taken literally; empty values are null, no forward-filling."""
        rows = [
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 1",
                "Concentration_Units": "microM",
                "MOA": "Inhibitor",
                "Data0": "10",
                "Conc0": "0.1",
            },
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 2",
                "Concentration_Units": "microM",
                "MOA": "",  # Empty; not filled from row 1
                "Data0": "20",
                "Conc0": "0.1",
            },
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, report = RawDataset.load_from_file(
                csv_file,
                skip_control_response_normalization=True,
                response_normalization="asymptote_normalized",
            )
            assert len(raw_data) == 2
            profiles = list(raw_data.profiles)
            assert (
                profiles[0].metadata is not None and profiles[0].metadata.get("MOA") == "Inhibitor"
            )
            # Row 2 has empty MOA -> null, not forward-filled
            assert (
                profiles[1].metadata is None
                or "MOA" not in (profiles[1].metadata or {})
                or profiles[1].metadata.get("MOA") in ("", None)
            )
        finally:
            csv_file.unlink()

    def test_metadata_extraction(self):
        rows = [
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 1",
                "Concentration_Units": "microM",
                "MOA": "Inhibitor",
                "drug targets": "Target1, Target2",
                "Data0": "10",
                "Conc0": "0.1",
            }
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, report = RawDataset.load_from_file(
                csv_file,
                skip_control_response_normalization=True,
                response_normalization="asymptote_normalized",
            )
            profile = list(raw_data.profiles)[0]
            assert profile.metadata is not None
            assert "MOA" in profile.metadata
            assert profile.metadata["MOA"] == "Inhibitor"
        finally:
            csv_file.unlink()

    def test_empty_file(self):
        csv_file = self.create_test_csv([])
        try:
            raw_data, report = RawDataset.load_from_file(
                csv_file,
                skip_control_response_normalization=True,
                response_normalization="asymptote_normalized",
            )
            assert len(raw_data) == 0
        finally:
            csv_file.unlink()

    def test_add_profile(self):
        assay = Assay(name="Test")
        raw_data = RawDataset(assay=assay, response_normalization="asymptote_normalized")
        compound = Compound(name="Test", drug_id="TEST")
        cell_line = CellLine(name="Test")
        profile = DoseResponseProfile(compound=compound, cell_line=cell_line, assay=assay)
        raw_data.add_profile(profile)
        assert len(raw_data) == 1

    def test_add_duplicate_profile(self):
        assay = Assay(name="Test")
        raw_data = RawDataset(assay=assay, response_normalization="asymptote_normalized")
        compound = Compound(name="Test", drug_id="TEST")
        cell_line = CellLine(name="Test")
        profile = DoseResponseProfile(compound=compound, cell_line=cell_line, assay=assay)
        raw_data.add_profile(profile)
        with pytest.raises(ValueError):
            raw_data.add_profile(profile)

    def test_to_screening_dataset(self):
        rows = [
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 1",
                "Concentration_Units": "microM",
                "Data0": "10",
                "Data1": "20",
                "Data2": "50",
                "Data3": "90",
                "Conc0": "0.1",
                "Conc1": "1",
                "Conc2": "10",
                "Conc3": "100",
            }
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, report = RawDataset.load_from_file(
                csv_file,
                skip_control_response_normalization=True,
                response_normalization="asymptote_normalized",
            )
            screening_data, process_report = raw_data.to_screening_dataset()
            assert isinstance(screening_data, ScreeningDataset)
            assert len(screening_data) == 1
            profile = list(screening_data.profiles)[0]
            assert profile.s_prime is not None
            assert profile.hill_params is not None
        finally:
            csv_file.unlink()

    def test_strict_load_response_scale_vs_asymptote_normalized(self, tmp_path):
        """DMSO path: import-time response_normalization selects pipeline before Hill fit."""
        csv_file = tmp_path / "strict.csv"
        with open(csv_file, "w", newline="") as f:
            w = csv.DictWriter(
                f,
                fieldnames=[
                    "Compound Name",
                    "Compound_ID",
                    "Cell_Line",
                    "Concentration_Units",
                    "Control_Response",
                    "Responses",
                    "Concentrations",
                ],
            )
            w.writeheader()
            w.writerow(
                {
                    "Compound Name": "Demo",
                    "Compound_ID": "D1",
                    "Cell_Line": "CL1",
                    "Concentration_Units": "microM",
                    "Control_Response": "35.3",
                    "Responses": "41.24,40.86,38.78,38.38,6.64,0.94,0,0",
                    "Concentrations": "0.003,0.01,0.03,0.1,0.3,1,3,10",
                }
            )

        raw_a, _ = RawDataset.load_from_file(
            csv_file,
            values_as="list",
            response_normalization="asymptote_normalized",
        )
        raw_b, _ = RawDataset.load_from_file(
            csv_file,
            values_as="list",
            response_normalization="response_scale",
        )
        sa, _ = raw_a.to_screening_dataset()
        sb, _ = raw_b.to_screening_dataset()
        pa = list(sa.profiles)[0]
        pb = list(sb.profiles)[0]
        assert pa.s_prime is not None and pb.s_prime is not None
        assert abs(pa.s_prime - pb.s_prime) > 0.05


class TestScreeningDataset:
    """Tests for ScreeningDataset."""

    @pytest.fixture
    def screening_dataset(self):
        assay = Assay(name="Test")
        dataset = ScreeningDataset(assay=assay)
        compound = Compound(name="Test", drug_id="TEST")
        cell_line1 = CellLine(name="Cell Line 1")
        cell_line2 = CellLine(name="Cell Line 2")

        # Create profiles with S' values
        profile1 = DoseResponseProfile(
            compound=compound,
            cell_line=cell_line1,
            assay=assay,
            concentrations=[0.1, 1, 10, 100],
            responses=[10, 20, 50, 90],
        )
        profile1.fit_and_calculate_s_prime()

        profile2 = DoseResponseProfile(
            compound=compound,
            cell_line=cell_line2,
            assay=assay,
            concentrations=[0.1, 1, 10, 100],
            responses=[5, 15, 60, 95],
        )
        profile2.fit_and_calculate_s_prime()

        dataset.add_profile(profile1)
        dataset.add_profile(profile2)
        return dataset

    def test_add_profile_without_s_prime(self):
        assay = Assay(name="Test")
        dataset = ScreeningDataset(assay=assay)
        compound = Compound(name="Test", drug_id="TEST")
        cell_line = CellLine(name="Test")
        profile = DoseResponseProfile(compound=compound, cell_line=cell_line, assay=assay)
        with pytest.raises(ValueError, match="must have S'"):
            dataset.add_profile(profile)

    def test_calculate_delta_s_prime_single(self, screening_dataset):
        results = screening_dataset.calculate_delta_s_prime(
            reference_cell_lines="Cell Line 1", test_cell_lines="Cell Line 2"
        )
        assert "Cell Line 1" in results
        assert len(results["Cell Line 1"]) == 1
        assert "delta_s_prime" in results["Cell Line 1"][0]

    def test_calculate_delta_s_prime_multiple(self, screening_dataset):
        results = screening_dataset.calculate_delta_s_prime(
            reference_cell_lines=["Cell Line 1"], test_cell_lines=["Cell Line 2"]
        )
        assert "Cell Line 1" in results

    def test_to_dict_list(self, screening_dataset):
        results = screening_dataset.to_dict_list()
        assert isinstance(results, list)
        assert len(results) == 2
        assert "compound_name" in results[0]
        assert "s_prime" in results[0]

    def test_export_to_csv(self, screening_dataset, tmp_path):
        output_file = tmp_path / "test_export.csv"
        screening_dataset.export_to_csv(output_file)
        assert output_file.exists()

        # Verify CSV content
        with open(output_file, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) == 2
            assert "Compound Name" in rows[0]
            assert "S'" in rows[0]


class TestSPrimeAPI:
    """Tests for SPrime main API class."""

    def test_load(self, tmp_path):
        csv_file = tmp_path / "test.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "Compound Name",
                    "Compound_ID",
                    "Cell_Line",
                    "Concentration_Units",
                    "Data0",
                    "Conc0",
                ],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "Compound Name": "Test",
                    "Compound_ID": "TEST",
                    "Cell_Line": "Cell1",
                    "Concentration_Units": "microM",
                    "Data0": "10",
                    "Conc0": "0.1",
                }
            )

        raw_data, report = SPrime.load(
            csv_file,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert isinstance(raw_data, RawDataset)
        assert len(raw_data) == 1
        assert report is not None

    def test_process(self, tmp_path):
        csv_file = tmp_path / "test.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "Compound Name",
                    "Compound_ID",
                    "Cell_Line",
                    "Concentration_Units",
                    "Data0",
                    "Data1",
                    "Data2",
                    "Data3",
                    "Conc0",
                    "Conc1",
                    "Conc2",
                    "Conc3",
                ],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "Compound Name": "Test",
                    "Compound_ID": "TEST",
                    "Cell_Line": "Cell1",
                    "Concentration_Units": "microM",
                    "Data0": "10",
                    "Data1": "30",
                    "Data2": "70",
                    "Data3": "90",
                    "Conc0": "0.1",
                    "Conc1": "1",
                    "Conc2": "10",
                    "Conc3": "100",
                }
            )

        raw_data, load_report = SPrime.load(
            csv_file,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        screening_data, process_report = SPrime.process(raw_data)
        assert isinstance(screening_data, ScreeningDataset)
        assert len(screening_data) == 1
        assert process_report is not None


class TestUtilityFunctions:
    """Tests for utility functions."""

    def test_fit_hill_from_raw_data(self):
        responses = [10, 20, 50, 90]
        concentrations = [0.1, 1, 10, 100]
        # Fully pre-scaled responses: no control ratio, no extra x100
        params = fit_hill_from_raw_data(
            responses,
            concentrations,
            skip_control_response_normalization=True,
            response_normalization="response_scale",
            scale_factor=1.0,
        )
        assert isinstance(params, HillCurveParams)
        assert params.ec50 > 0

    def test_fit_hill_from_raw_data_requires_control_when_not_skipping(self):
        with pytest.raises(ValueError, match="control_response is required"):
            fit_hill_from_raw_data(
                [41.24, 40.86],
                [0.003, 0.01],
                skip_control_response_normalization=False,
                response_normalization="response_scale",
            )
        with pytest.raises(ValueError, match="non-zero"):
            fit_hill_from_raw_data(
                [41.24, 40.86],
                [0.003, 0.01],
                control_response=0.0,
                skip_control_response_normalization=False,
                response_normalization="response_scale",
            )

    def test_calculate_s_prime_from_params(self):
        # Inhibitor-style ordering: high response at x->0, low at saturation => Zero > Inf => S' > 0
        s_prime = calculate_s_prime_from_params(ac50=10.0, zero_asymptote=100.0, inf_asymptote=0.0)
        assert isinstance(s_prime, float)
        assert s_prime > 0

    def test_specific_s_prime_value_asymptote_normalized(self):
        """
        Golden: variation reference **Parameter Values** (normalized path) publish **A-D (span)** and
        **C (EC50)** only. S' = asinh(span/EC50); do **not** assume ``inf_asymptote`` (e.g. 0) - use
        raw->fit for actual asymptotes (see ``test_s_prime_from_raw_data_asymptote_normalized``).
        """
        ec50 = 0.202
        span = 97.0276  # A-D (Emax) = zero_asymptote - inf_asymptote for that fit
        s_prime = math.asinh(span / ec50)
        assert s_prime == pytest.approx(6.867631319810275, rel=1e-12)
        assert round(s_prime, 2) == 6.87

    def test_specific_s_prime_value_response_scale(self):
        """
        Variation reference **Parameter Values** (response-scale / x100 path): tabulated **A-D** is
        the span ``zero_asymptote - inf_asymptote``; **C** is EC50. For the published rounded row,
        S' = asinh(span/EC50). Actual ``zero_asymptote`` / ``inf_asymptote`` come only from **fit**
        or **precalc columns** - never assume (e.g.) ``inf_asymptote = 0``; see
        ``test_s_prime_from_raw_data_response_scale``.
        """
        ec50 = 0.20
        span = 113.34  # A-D (Emax) from the sheet for this path
        s_prime = math.asinh(span / ec50)
        assert s_prime == pytest.approx(7.032978022189737, rel=1e-9)
        assert (span / ec50) == pytest.approx(566.7, rel=1e-12)

    def test_s_prime_from_raw_data_asymptote_normalized(self):
        """
        Golden test: S' from raw data using asymptote-normalized process.
        Raw data from SPrime_variation_reference.csv (17-AAG, NF2-/-, 48h): concentrations and
        percent nucleus response; DMSO control 35.3. Process: (1) DMSO-relative ratio,
        (2) normalize_to_max_value on ratios (max -> 1), (3) x100 scale,
        (4) fit 4-PL, (5) S' from fitted params. Expect S' rounds to 6.87.
        """
        # Raw data from reference CSV (rows 2-9, exclude DMSO row)
        concentrations = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]  # microM
        raw_responses = [41.24, 40.86, 38.78, 38.38, 6.64, 0.94, 0.00, 0.00]  # percent nucleus
        dmso_control = 35.3

        # Integrated API (DMSO + asymptote_normalized) matches pipeline + legacy fit.
        params = fit_hill_from_raw_data(
            raw_responses,
            concentrations,
            control_response=dmso_control,
            skip_control_response_normalization=False,
            response_normalization="asymptote_normalized",
        )
        legacy = fit_hill_from_raw_data(
            pipeline_asymptote_normalized(raw_responses, dmso_control),
            concentrations,
            skip_control_response_normalization=True,
            response_normalization="response_scale",
            scale_factor=1.0,
        )
        assert params.ec50 == pytest.approx(legacy.ec50, rel=1e-9)
        assert params.zero_asymptote == pytest.approx(legacy.zero_asymptote, rel=1e-9)
        assert params.inf_asymptote == pytest.approx(legacy.inf_asymptote, rel=1e-9)

        s_prime = calculate_s_prime_from_params(
            params.ec50, params.zero_asymptote, params.inf_asymptote
        )
        assert s_prime == pytest.approx(6.867631094104203, rel=1e-9)

    def test_s_prime_from_raw_data_response_scale(self):
        """
        Golden test: S' from raw data using response-scale (non-normalized) process.
        Same raw data and DMSO control. Process: (1) DMSO-relative ratio, (2) x100, no max normalization,
        (3) fit 4-PL, (4) S' from fitted params. Tabulated reference **A-D** is the span
        **(zero_asymptote - inf_asymptote)** with **inf_asymptote ~= 0.85** (not 0) on the response
        scale; see ``test_specific_s_prime_value_response_scale``. The spreadsheet **Parameter Values**
        row gives a slightly different S' (**~7.033**) from rounded EC50/span; **auto-fit** here is
        anchored at **7.022989767** (``abs=2e-4`` allows tiny ``curve_fit`` / platform drift).
        """
        concentrations = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
        raw_responses = [41.24, 40.86, 38.78, 38.38, 6.64, 0.94, 0.00, 0.00]
        dmso_control = 35.3

        params = fit_hill_from_raw_data(
            raw_responses,
            concentrations,
            control_response=dmso_control,
            skip_control_response_normalization=False,
            response_normalization="response_scale",
        )
        legacy = fit_hill_from_raw_data(
            pipeline_response_scale(raw_responses, dmso_control),
            concentrations,
            skip_control_response_normalization=True,
            response_normalization="response_scale",
            scale_factor=1.0,
        )
        assert params.ec50 == pytest.approx(legacy.ec50, rel=1e-9)

        s_prime = calculate_s_prime_from_params(
            params.ec50, params.zero_asymptote, params.inf_asymptote
        )
        # Nominal non-normalized auto-fit S' (reference checkpoint); see docstring.
        assert s_prime == pytest.approx(7.022989767, abs=2e-4)

    def test_s_prime_from_precomputed_ratios_matches_dmso_paths(self):
        """Non-DMSO-in-API: pass test/control ratios from the same slice as variation reference."""
        concentrations = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
        raw_responses = [41.24, 40.86, 38.78, 38.38, 6.64, 0.94, 0.00, 0.00]
        dmso = 35.3
        ratios = [r / dmso for r in raw_responses]

        p_dmso_norm = fit_hill_from_raw_data(
            raw_responses,
            concentrations,
            control_response=dmso,
            skip_control_response_normalization=False,
            response_normalization="asymptote_normalized",
        )
        p_ratio_norm = fit_hill_from_raw_data(
            ratios,
            concentrations,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert p_dmso_norm.ec50 == pytest.approx(p_ratio_norm.ec50, rel=1e-9)
        assert calculate_s_prime_from_params(
            p_dmso_norm.ec50, p_dmso_norm.zero_asymptote, p_dmso_norm.inf_asymptote
        ) == pytest.approx(
            calculate_s_prime_from_params(
                p_ratio_norm.ec50, p_ratio_norm.zero_asymptote, p_ratio_norm.inf_asymptote
            ),
            rel=1e-9,
        )

        p_dmso_rs = fit_hill_from_raw_data(
            raw_responses,
            concentrations,
            control_response=dmso,
            skip_control_response_normalization=False,
            response_normalization="response_scale",
        )
        p_ratio_rs = fit_hill_from_raw_data(
            ratios,
            concentrations,
            skip_control_response_normalization=True,
            response_normalization="response_scale",
        )
        assert p_dmso_rs.ec50 == pytest.approx(p_ratio_rs.ec50, rel=1e-9)
        assert calculate_s_prime_from_params(
            p_dmso_rs.ec50, p_dmso_rs.zero_asymptote, p_dmso_rs.inf_asymptote
        ) == pytest.approx(
            calculate_s_prime_from_params(
                p_ratio_rs.ec50, p_ratio_rs.zero_asymptote, p_ratio_rs.inf_asymptote
            ),
            rel=1e-9,
        )

    def test_fit_hill_from_raw_data_concentration_units_nM_matches_uM(self):
        raw_responses = [41.24, 40.86, 38.78, 38.38, 6.64, 0.94, 0.00, 0.00]
        conc_uM = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]
        conc_nM = [c * 1000.0 for c in conc_uM]
        dmso = 35.3
        p_um = fit_hill_from_raw_data(
            raw_responses,
            conc_uM,
            control_response=dmso,
            skip_control_response_normalization=False,
            response_normalization="response_scale",
        )
        p_nm = fit_hill_from_raw_data(
            raw_responses,
            conc_nM,
            control_response=dmso,
            skip_control_response_normalization=False,
            response_normalization="response_scale",
            concentration_units="nM",
        )
        assert p_um.ec50 == pytest.approx(p_nm.ec50, rel=1e-9)

    def test_get_s_primes_from_file(self, tmp_path):
        csv_file = tmp_path / "test.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "Compound Name",
                    "Compound_ID",
                    "Cell_Line",
                    "Concentration_Units",
                    "Data0",
                    "Data1",
                    "Data2",
                    "Data3",
                    "Conc0",
                    "Conc1",
                    "Conc2",
                    "Conc3",
                ],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "Compound Name": "Test",
                    "Compound_ID": "TEST",
                    "Cell_Line": "Cell1",
                    "Concentration_Units": "microM",
                    "Data0": "10",
                    "Data1": "30",
                    "Data2": "70",
                    "Data3": "90",
                    "Conc0": "0.1",
                    "Conc1": "1",
                    "Conc2": "10",
                    "Conc3": "100",
                }
            )

        results = get_s_primes_from_file(
            csv_file,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert isinstance(results, list)
        assert len(results) == 1
        assert "s_prime" in results[0]

    def test_get_s_prime_from_data(self):
        list_of_rows = [
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 1",
                "Concentration_Units": "microM",
                "Data0": "10",
                "Data1": "20",
                "Data2": "50",
                "Data3": "90",
                "Conc0": "0.1",
                "Conc1": "1",
                "Conc2": "10",
                "Conc3": "100",
            }
        ]
        results = get_s_prime_from_data(
            list_of_rows,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert isinstance(results, list)
        assert len(results) == 1
        assert "s_prime" in results[0]

    def test_get_s_prime_from_data_values_as_list(self):
        """get_s_prime_from_data with values_as='list' and Responses/Concentrations."""
        list_of_rows = [
            {
                "Compound Name": "Drug A",
                "Compound_ID": "DRUG001",
                "Cell_Line": "Cell Line 1",
                "Concentration_Units": "microM",
                "Responses": "10,20,50,90",
                "Concentrations": "0.1,1,10,100",
            }
        ]
        results = get_s_prime_from_data(
            list_of_rows,
            values_as="list",
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert isinstance(results, list)
        assert len(results) == 1
        assert "s_prime" in results[0]
        assert results[0]["s_prime"] is not None

    def test_calculate_delta_s_prime_with_screening_dataset(self):
        # Create a simple screening dataset
        assay = Assay(name="Test")
        dataset = ScreeningDataset(assay=assay)
        compound = Compound(name="Test", drug_id="TEST")

        for i, cell_line_name in enumerate(["Ref", "Test"]):
            cell_line = CellLine(name=cell_line_name)
            profile = DoseResponseProfile(
                compound=compound,
                cell_line=cell_line,
                assay=assay,
                concentrations=[0.1, 1, 10, 100],
                responses=[10 + i * 10, 20 + i * 10, 50 + i * 10, 90 + i * 10],
            )
            profile.fit_and_calculate_s_prime()
            dataset.add_profile(profile)

        results = calculate_delta_s_prime(
            dataset, reference_cell_line_names="Ref", test_cell_line_names="Test"
        )
        assert "Ref" in results
        assert len(results["Ref"]) == 1

    def test_calculate_delta_s_prime_with_list_dict(self):
        # Create list of dicts with S' values
        list_of_dicts = [
            {
                "compound_name": "Test",
                "drug_id": "TEST",
                "cell_line": "Ref",
                "s_prime": 2.5,
                "ec50": 10.0,
                "zero_asymptote": 0.0,
                "inf_asymptote": 100.0,
                "steepness_coefficient": 1.5,
                "r_squared": 0.95,
            },
            {
                "compound_name": "Test",
                "drug_id": "TEST",
                "cell_line": "Test",
                "s_prime": 3.0,
                "ec50": 8.0,
                "zero_asymptote": 0.0,
                "inf_asymptote": 100.0,
                "steepness_coefficient": 1.5,
                "r_squared": 0.95,
            },
        ]

        results = calculate_delta_s_prime(
            list_of_dicts, reference_cell_line_names="Ref", test_cell_line_names="Test"
        )
        assert "Ref" in results
        assert len(results["Ref"]) == 1
        assert results["Ref"][0]["delta_s_prime"] == -0.5  # 2.5 - 3.0

    def test_convert_to_micromolar(self):
        # Test various unit conversions (smallest to largest)
        assert convert_to_micromolar([1.0], "fM") == [1e-9]
        assert convert_to_micromolar([1.0], "pM") == [1e-6]
        assert convert_to_micromolar([1.0], "nM") == [0.001]
        assert convert_to_micromolar([1.0], "microM") == [1.0]
        assert convert_to_micromolar([1.0], "mM") == [1000.0]
        assert convert_to_micromolar([1.0], "M") == [1000000.0]


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_csv_missing_required_columns(self, tmp_path):
        csv_file = tmp_path / "test.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["Some Column"])
            writer.writeheader()
            writer.writerow({"Some Column": "value"})

        # Should raise ValueError for missing required columns
        with pytest.raises(ValueError, match="Required column 'Cell_Line' not found"):
            raw_data, report = RawDataset.load_from_file(
                csv_file,
                skip_control_response_normalization=True,
                response_normalization="asymptote_normalized",
            )

    def test_csv_inconsistent_column_naming(self, tmp_path):
        # Test both Data0 and DATA0 formats
        csv_file = tmp_path / "test.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "Compound Name",
                    "Compound_ID",
                    "Cell_Line",
                    "Concentration_Units",
                    "DATA0",
                    "CONC0",
                ],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "Compound Name": "Test",
                    "Compound_ID": "TEST",
                    "Cell_Line": "Cell1",
                    "Concentration_Units": "microM",
                    "DATA0": "10",
                    "CONC0": "0.1",
                }
            )

        raw_data, report = RawDataset.load_from_file(
            csv_file,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert len(raw_data) == 1

    def test_column_order_does_not_matter(self, tmp_path):
        """Column order in CSV does not matter; matching is by heading name only."""
        row = {
            "Compound Name": "Drug X",
            "Compound_ID": "DX",
            "Cell_Line": "CL1",
            "Concentration_Units": "microM",
            "Data0": "5",
            "Data1": "15",
            "Data2": "40",
            "Data3": "85",
            "Conc0": "0.01",
            "Conc1": "0.1",
            "Conc2": "1",
            "Conc3": "10",
        }
        order_a = [
            "Compound_ID",
            "Cell_Line",
            "Concentration_Units",
            "Data0",
            "Data1",
            "Data2",
            "Data3",
            "Conc0",
            "Conc1",
            "Conc2",
            "Conc3",
            "Compound Name",
        ]
        order_b = [
            "Compound Name",
            "Data3",
            "Data2",
            "Data1",
            "Data0",
            "Conc3",
            "Conc2",
            "Conc1",
            "Conc0",
            "Cell_Line",
            "Concentration_Units",
            "Compound_ID",
        ]

        def write_csv(path, fieldnames):
            with open(path, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
                w.writeheader()
                w.writerow(row)

        path_a = tmp_path / "order_a.csv"
        path_b = tmp_path / "order_b.csv"
        write_csv(path_a, order_a)
        write_csv(path_b, order_b)

        raw_a, _ = RawDataset.load_from_file(
            path_a,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        raw_b, _ = RawDataset.load_from_file(
            path_b,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert len(raw_a) == 1 and len(raw_b) == 1
        p_a = list(raw_a.profiles)[0]
        p_b = list(raw_b.profiles)[0]
        assert p_a.compound.drug_id == p_b.compound.drug_id == "DX"
        assert p_a.cell_line.name == p_b.cell_line.name == "CL1"
        assert p_a.compound.name == p_b.compound.name == "Drug X"
        assert p_a.concentrations == p_b.concentrations
        assert p_a.responses == p_b.responses

        from sprime import SPrime as sp

        sa, _ = sp.process(raw_a)
        sb, _ = sp.process(raw_b)
        pa = list(sa.profiles)[0]
        pb = list(sb.profiles)[0]
        assert pa.s_prime is not None and pb.s_prime is not None
        assert abs(pa.s_prime - pb.s_prime) < 1e-9

    def test_csv_scientific_notation(self, tmp_path):
        csv_file = tmp_path / "test.csv"
        with open(csv_file, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "Compound Name",
                    "Compound_ID",
                    "Cell_Line",
                    "Concentration_Units",
                    "Data0",
                    "Conc0",
                ],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "Compound Name": "Test",
                    "Compound_ID": "TEST",
                    "Cell_Line": "Cell1",
                    "Concentration_Units": "microM",
                    "Data0": "10",
                    "Conc0": "1.30E-09",  # Scientific notation
                }
            )

        raw_data, report = RawDataset.load_from_file(
            csv_file,
            skip_control_response_normalization=True,
            response_normalization="asymptote_normalized",
        )
        assert len(raw_data) == 1
        profile = list(raw_data.profiles)[0]
        assert len(profile.concentrations) == 1

    def test_insufficient_data_points(self):
        compound = Compound(name="Test", drug_id="TEST")
        cell_line = CellLine(name="Test")
        assay = Assay(name="Test")
        profile = DoseResponseProfile(
            compound=compound,
            cell_line=cell_line,
            assay=assay,
            concentrations=[0.1, 1, 10],  # Only 3 points
            responses=[10, 20, 30],
        )
        with pytest.raises((ValueError, RuntimeError)):
            profile.fit_hill_curve()

    def test_delta_s_prime_missing_profiles(self):
        assay = Assay(name="Test")
        dataset = ScreeningDataset(assay=assay)
        compound = Compound(name="Test", drug_id="TEST")
        cell_line = CellLine(name="Ref")

        profile = DoseResponseProfile(
            compound=compound,
            cell_line=cell_line,
            assay=assay,
            concentrations=[0.1, 1, 10, 100],
            responses=[10, 20, 50, 90],
        )
        profile.fit_and_calculate_s_prime()
        dataset.add_profile(profile)

        # Try to calculate delta S' with missing test cell line
        results = dataset.calculate_delta_s_prime(
            reference_cell_lines="Ref",
            test_cell_lines="Test",  # Doesn't exist
        )
        # Should return empty results, not crash
        assert "Ref" in results
        assert len(results["Ref"]) == 0
