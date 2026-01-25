"""
Integration tests for sprime.

Tests complete workflows end-to-end using real data examples.
"""

import pytest
import sys
import tempfile
import csv
from pathlib import Path

from sprime import (
    SPrime, RawDataset, ScreeningDataset,
    get_s_primes_from_file, calculate_delta_s_prime
)


class TestCompleteWorkflow:
    """Test complete workflow: Load → Process → Calculate S' → Delta S'."""
    
    def create_test_csv(self, rows):
        """Helper to create temporary CSV file."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv', newline='')
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())
        
        writer = csv.DictWriter(temp_file, fieldnames=sorted(fieldnames))
        writer.writeheader()
        writer.writerows(rows)
        temp_file.close()
        return Path(temp_file.name)
    
    def test_complete_workflow_single_compound_multiple_cell_lines(self):
        """Test complete workflow with one compound tested against multiple cell lines."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Reference Cell Line',
                'Concentration_Units': 'microM',
                'Data0': '5',
                'Data1': '10',
                'Data2': '30',
                'Data3': '70',
                'Data4': '85',
                'Conc0': '0.01',
                'Conc1': '0.1',
                'Conc2': '1',
                'Conc3': '10',
                'Conc4': '100',
            },
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Test Cell Line',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Data4': '95',
                'Conc0': '0.01',
                'Conc1': '0.1',
                'Conc2': '1',
                'Conc3': '10',
                'Conc4': '100',
            }
        ]
        
        csv_file = self.create_test_csv(rows)
        try:
            # Step 1: Load
            raw_data, load_report = SPrime.load(csv_file)
            assert len(raw_data) == 2
            
            # Step 2: Process
            screening_data, process_report = SPrime.process(raw_data)
            assert len(screening_data) == 2
            
            # Verify all profiles have S' calculated
            for profile in screening_data.profiles:
                assert profile.s_prime is not None
                assert profile.hill_params is not None
                assert profile.hill_params.r_squared is not None
            
            # Step 3: Calculate Delta S'
            delta_results = screening_data.calculate_delta_s_prime(
                reference_cell_lines="Reference Cell Line",
                test_cell_lines="Test Cell Line"
            )
            
            assert "Reference Cell Line" in delta_results
            assert len(delta_results["Reference Cell Line"]) == 1
            assert "delta_s_prime" in delta_results["Reference Cell Line"][0]
            
            # Step 4: Export
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
                export_file = Path(f.name)
            
            screening_data.export_to_csv(export_file)
            assert export_file.exists()
            
            # Verify export content
            with open(export_file, 'r') as f:
                reader = csv.DictReader(f)
                exported_rows = list(reader)
                assert len(exported_rows) == 2
            
            export_file.unlink()
            
        finally:
            csv_file.unlink()
    
    def test_complete_workflow_multiple_compounds(self):
        """Test workflow with multiple compounds."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Cell Line 1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            },
            {
                'Compound Name': 'Drug B',
                'Compound_ID': 'DRUG002',
                'Cell_Line': 'Cell Line 1',
                'Concentration_Units': 'microM',
                'Data0': '5',
                'Data1': '15',
                'Data2': '60',
                'Data3': '95',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            }
        ]
        
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, load_report = SPrime.load(csv_file)
            screening_data, process_report = SPrime.process(raw_data)
            
            assert len(screening_data) == 2
            
            # Verify we have results for both compounds
            results = screening_data.to_dict_list()
            assert len(results) == 2
            # Verify both have S' values calculated
            s_primes = [r['s_prime'] for r in results if r['s_prime'] is not None]
            assert len(s_primes) == 2
            # Verify both compounds are present
            drug_ids = {r['drug_id'] for r in results}
            assert 'DRUG001' in drug_ids
            assert 'DRUG002' in drug_ids
            
        finally:
            csv_file.unlink()
    
    def test_workflow_with_precalculated_params(self):
        """Test workflow with pre-calculated Hill parameters."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Cell Line 1',
                'AC50': '10.0',
                'Upper': '100.0',
                'Lower': '0.0',
                'Hill_Slope': '1.5',
                'r2': '0.95',
            }
        ]
        
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, load_report = SPrime.load(csv_file)
            assert len(raw_data) == 1
            
            profile = list(raw_data.profiles)[0]
            assert profile.hill_params is not None
            assert profile.hill_params.ec50 == 10.0
            
            screening_data, process_report = SPrime.process(raw_data)
            assert len(screening_data) == 1
            
            profile = list(screening_data.profiles)[0]
            assert profile.s_prime is not None
            
        finally:
            csv_file.unlink()
    
    def test_workflow_with_metadata(self):
        """Test workflow preserving metadata."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Cell Line 1',
                'Concentration_Units': 'microM',
                'MOA': 'Inhibitor',
                'drug targets': 'Target1, Target2',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            }
        ]
        
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, load_report = SPrime.load(csv_file)
            screening_data, process_report = SPrime.process(raw_data)
            
            profile = list(screening_data.profiles)[0]
            assert profile.metadata is not None
            assert 'MOA' in profile.metadata
            assert profile.metadata['MOA'] == 'Inhibitor'
            
            # Export and verify metadata is included
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
                export_file = Path(f.name)
            
            screening_data.export_to_csv(export_file, include_metadata=True)
            
            with open(export_file, 'r') as f:
                reader = csv.DictReader(f)
                exported_rows = list(reader)
                assert 'MOA' in exported_rows[0]
                assert exported_rows[0]['MOA'] == 'Inhibitor'
            
            export_file.unlink()
            
        finally:
            csv_file.unlink()
    
    def test_delta_s_prime_export(self):
        """Test delta S' calculation and export."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Reference',
                'Concentration_Units': 'microM',
                'Data0': '5',
                'Data1': '10',
                'Data2': '30',
                'Data3': '70',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            },
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Test',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            }
        ]
        
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, load_report = SPrime.load(csv_file)
            screening_data, process_report = SPrime.process(raw_data)
            
            delta_results = screening_data.calculate_delta_s_prime(
                reference_cell_lines="Reference",
                test_cell_lines="Test"
            )
            
            # Export delta S' results
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
                export_file = Path(f.name)
            
            ScreeningDataset.export_delta_s_prime_to_csv(delta_results, export_file)
            assert export_file.exists()
            
            # Verify export content
            with open(export_file, 'r') as f:
                reader = csv.DictReader(f)
                exported_rows = list(reader)
                assert len(exported_rows) == 1
                assert "Delta S'" in exported_rows[0]
                assert "Rank" in exported_rows[0]
            
            export_file.unlink()
            
        finally:
            csv_file.unlink()

    def test_assay_timescale_propagates_raw(self):
        """Raw import (DATA/CONC): assay_timescale 24hr/48hr propagates to S' metadata and master export."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Reference',
                'Concentration_Units': 'microM',
                'assay_timescale': '24hr',
                'Data0': '5',
                'Data1': '10',
                'Data2': '30',
                'Data3': '70',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            },
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Test',
                'Concentration_Units': 'microM',
                'assay_timescale': '48hr',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            },
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, _ = SPrime.load(csv_file)
            screening_data, _ = SPrime.process(raw_data)
            profiles = sorted(screening_data.profiles, key=lambda p: p.cell_line.name)
            assert len(profiles) == 2
            ref_prof = next(p for p in profiles if p.cell_line.name == 'Reference')
            test_prof = next(p for p in profiles if p.cell_line.name == 'Test')
            assert ref_prof.metadata is not None and ref_prof.metadata.get('assay_timescale') == '24hr'
            assert test_prof.metadata is not None and test_prof.metadata.get('assay_timescale') == '48hr'

            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
                export_file = Path(f.name)
            screening_data.export_to_csv(export_file, include_metadata=True)
            with open(export_file, 'r') as f:
                reader = csv.DictReader(f)
                exported = list(reader)
            assert 'assay_timescale' in exported[0]
            vals = {r['Cell_Line']: r['assay_timescale'] for r in exported}
            assert vals['Reference'] == '24hr' and vals['Test'] == '48hr'
            export_file.unlink()
        finally:
            csv_file.unlink()

    def test_assay_timescale_propagates_precalc(self):
        """Precalc import (AC50/Upper/Lower only): assay_timescale 24hr/48hr propagates to S' metadata and master export."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Reference',
                'AC50': '10.0',
                'Upper': '100.0',
                'Lower': '0.0',
                'assay_timescale': '24hr',
            },
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Test',
                'AC50': '8.0',
                'Upper': '100.0',
                'Lower': '0.0',
                'assay_timescale': '48hr',
            },
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, _ = SPrime.load(csv_file)
            screening_data, _ = SPrime.process(raw_data)
            profiles = sorted(screening_data.profiles, key=lambda p: p.cell_line.name)
            assert len(profiles) == 2
            ref_prof = next(p for p in profiles if p.cell_line.name == 'Reference')
            test_prof = next(p for p in profiles if p.cell_line.name == 'Test')
            assert ref_prof.metadata is not None and ref_prof.metadata.get('assay_timescale') == '24hr'
            assert test_prof.metadata is not None and test_prof.metadata.get('assay_timescale') == '48hr'

            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
                export_file = Path(f.name)
            screening_data.export_to_csv(export_file, include_metadata=True)
            with open(export_file, 'r') as f:
                reader = csv.DictReader(f)
                exported = list(reader)
            assert 'assay_timescale' in exported[0]
            vals = {r['Cell_Line']: r['assay_timescale'] for r in exported}
            assert vals['Reference'] == '24hr' and vals['Test'] == '48hr'
            export_file.unlink()
        finally:
            csv_file.unlink()

    def test_delta_s_prime_output_assay_timescale(self):
        """Delta S' with headings_one_to_one_in_ref_and_test=['assay_timescale']; delta export includes assay_timescale."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Reference',
                'Concentration_Units': 'microM',
                'assay_timescale': '24hr',
                'Data0': '5',
                'Data1': '10',
                'Data2': '30',
                'Data3': '70',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            },
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Test',
                'Concentration_Units': 'microM',
                'assay_timescale': '48hr',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            },
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, _ = SPrime.load(csv_file)
            screening_data, _ = SPrime.process(raw_data)
            delta_results = screening_data.calculate_delta_s_prime(
                reference_cell_lines='Reference',
                test_cell_lines='Test',
                headings_one_to_one_in_ref_and_test=['assay_timescale'],
                source_profile='test',
            )
            assert 'Reference' in delta_results
            assert len(delta_results['Reference']) == 1
            comp = delta_results['Reference'][0]
            assert 'assay_timescale' in comp
            assert comp['assay_timescale'] == '48hr'

            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as f:
                export_file = Path(f.name)
            ScreeningDataset.export_delta_s_prime_to_csv(
                delta_results,
                export_file,
                headings_one_to_one_in_ref_and_test=['assay_timescale'],
            )
            with open(export_file, 'r') as f:
                reader = csv.DictReader(f)
                exported = list(reader)
            assert len(exported) == 1
            assert 'assay_timescale' in exported[0]
            assert exported[0]['assay_timescale'] == '48hr'
            export_file.unlink()
        finally:
            csv_file.unlink()


class TestAllowOverwriteHillCoefficients:
    """Test allow_overwrite_hill_coefficients: raise when disallowed, succeed + warn when allowed."""

    def create_test_csv(self, rows):
        """Helper to create temporary CSV file."""
        temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv', newline='')
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())
        writer = csv.DictWriter(temp_file, fieldnames=sorted(fieldnames))
        writer.writeheader()
        writer.writerows(rows)
        temp_file.close()
        return Path(temp_file.name)

    def test_allow_overwrite_hill_coefficients_raises(self):
        """With raw + pre-calc, process(allow_overwrite=False) raises ValueError."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Cell Line 1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
                'AC50': '10.0',
                'Upper': '100.0',
                'Lower': '0.0',
                'Hill_Slope': '1.5',
                'r2': '0.95',
            }
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, _ = SPrime.load(csv_file)
            assert len(raw_data) == 1
            with pytest.raises(ValueError, match="allow_overwrite_hill_coefficients"):
                SPrime.process(raw_data, allow_overwrite_hill_coefficients=False)
        finally:
            csv_file.unlink()

    def test_allow_overwrite_hill_coefficients_succeeds(self):
        """With raw + pre-calc, process(allow_overwrite=True) succeeds and logs OVERWRITE_HILL."""
        rows = [
            {
                'Compound Name': 'Drug A',
                'Compound_ID': 'DRUG001',
                'Cell_Line': 'Cell Line 1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
                'AC50': '10.0',
                'Upper': '100.0',
                'Lower': '0.0',
                'Hill_Slope': '1.5',
                'r2': '0.95',
            }
        ]
        csv_file = self.create_test_csv(rows)
        try:
            raw_data, _ = SPrime.load(csv_file)
            screening_data, process_report = SPrime.process(
                raw_data, allow_overwrite_hill_coefficients=True
            )
            assert len(screening_data) == 1
            profile = list(screening_data.profiles)[0]
            assert profile.hill_params is not None
            assert profile.s_prime is not None
            overwrite_warnings = [w for w in process_report.warnings if w.category == "OVERWRITE_HILL"]
            assert len(overwrite_warnings) >= 1
        finally:
            csv_file.unlink()


class TestErrorHandling:
    """Test error handling in workflows."""
    
    def test_workflow_with_failed_curve_fit(self, tmp_path):
        """Test handling of profiles that fail curve fitting."""
        csv_file = tmp_path / "test.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Compound Name', 'Compound_ID', 'Cell_Line', 'Concentration_Units', 'Data0', 'Conc0'])
            writer.writeheader()
            writer.writerow({
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Conc0': '0.1'
                # Only 1 data point - will fail fitting
            })
        
        raw_data, load_report = SPrime.load(csv_file)
        
        # Processing should raise ValueError when curve fitting fails
        with pytest.raises(ValueError, match="Curve fitting failed"):
            screening_data, process_report = SPrime.process(raw_data)
    
    def test_workflow_with_missing_cell_line(self, tmp_path):
        """Test delta S' calculation with missing cell line."""
        csv_file = tmp_path / "test.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Compound Name', 'Compound_ID', 'Cell_Line', 'Concentration_Units', 'Data0', 'Data1', 'Data2', 'Data3', 'Conc0', 'Conc1', 'Conc2', 'Conc3'])
            writer.writeheader()
            writer.writerow({
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '30',
                'Data2': '70',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100'
            })
        
        raw_data, load_report = SPrime.load(csv_file)
        screening_data, process_report = SPrime.process(raw_data)
        
        # Try to calculate delta S' with non-existent cell line
        results = screening_data.calculate_delta_s_prime(
            reference_cell_lines="Cell1",
            test_cell_lines="NonExistent"
        )
        
        # Should return empty results, not crash
        assert "Cell1" in results
        assert len(results["Cell1"]) == 0

    def test_raw_without_concentration_units_raises(self, tmp_path):
        """Raw CSV with DATA/CONC but no Concentration_Units column raises ValueError."""
        csv_file = tmp_path / "test.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Compound Name', 'Compound_ID', 'Cell_Line', 'Data0', 'Data1', 'Data2', 'Data3', 'Conc0', 'Conc1', 'Conc2', 'Conc3'])
            writer.writeheader()
            writer.writerow({
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Data0': '10',
                'Data1': '30',
                'Data2': '70',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100',
            })
        with pytest.raises(ValueError, match="Required column 'Concentration_Units' not found"):
            SPrime.load(csv_file)

    def test_values_as_list_load_process(self):
        """Load demo_data_raw_list with values_as='list', process, check S'."""
        from pathlib import Path
        csv_path = Path(__file__).resolve().parent.parent / "docs" / "usage" / "demo_data_raw_list.csv"
        raw, _ = SPrime.load(csv_path, values_as="list")
        assert len(raw) == 1
        screening, _ = SPrime.process(raw)
        assert len(screening) == 1
        p = list(screening.profiles)[0]
        assert p.s_prime is not None
        assert p.hill_params is not None

    def test_values_as_list_missing_responses_raises(self, tmp_path):
        """CSV with Concentrations but no Responses, values_as='list' -> ValueError."""
        csv_file = tmp_path / "test.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(
                f,
                fieldnames=['Compound Name', 'Compound_ID', 'Cell_Line', 'Concentration_Units', 'Concentrations']
            )
            writer.writeheader()
            writer.writerow({
                'Compound Name': 'Test',
                'Compound_ID': 'T',
                'Cell_Line': 'C',
                'Concentration_Units': 'microM',
                'Concentrations': '1,2,3,4',
            })
        with pytest.raises(ValueError, match="No dose-response data"):
            SPrime.load(csv_file, values_as="list")

    def test_values_as_list_mismatched_lengths_raises(self, tmp_path):
        """Responses and Concentrations length mismatch -> ValueError."""
        csv_file = tmp_path / "test.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(
                f,
                fieldnames=['Compound Name', 'Compound_ID', 'Cell_Line', 'Concentration_Units', 'Responses', 'Concentrations']
            )
            writer.writeheader()
            writer.writerow({
                'Compound Name': 'Test',
                'Compound_ID': 'T',
                'Cell_Line': 'C',
                'Concentration_Units': 'microM',
                'Responses': '1,2,3,4',
                'Concentrations': '1,2,3',
            })
        with pytest.raises(ValueError, match="length mismatch"):
            SPrime.load(csv_file, values_as="list")


class TestGetSPrimesFromFile:
    """Test get_s_primes_from_file utility function."""
    
    def test_get_s_primes_from_file_complete(self, tmp_path):
        """Test get_s_primes_from_file with complete data."""
        csv_file = tmp_path / "test.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['Compound Name', 'Compound_ID', 'Cell_Line', 'Concentration_Units', 'Data0', 'Data1', 'Data2', 'Data3', 'Conc0', 'Conc1', 'Conc2', 'Conc3'])
            writer.writeheader()
            writer.writerow({
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '30',
                'Data2': '70',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100'
            })
        
        results = get_s_primes_from_file(csv_file)
        assert isinstance(results, list)
        assert len(results) == 1
        assert "s_prime" in results[0]
        assert "rank" in results[0]
        assert results[0]["rank"] == 1


class TestCalculateDeltaSPrimeUtility:
    """Test calculate_delta_s_prime utility function."""
    
    def test_calculate_delta_s_prime_with_list_dict(self):
        """Test calculate_delta_s_prime with list of dicts."""
        list_of_dicts = [
            {
                'compound_name': 'Drug A',
                'drug_id': 'DRUG001',
                'cell_line': 'Reference',
                's_prime': 2.0,
                'ec50': 10.0,
                'upper': 100.0,
                'lower': 0.0,
                'hill_coefficient': 1.5,
                'r_squared': 0.95
            },
            {
                'compound_name': 'Drug A',
                'drug_id': 'DRUG001',
                'cell_line': 'Test',
                's_prime': 3.0,
                'ec50': 8.0,
                'upper': 100.0,
                'lower': 0.0,
                'hill_coefficient': 1.5,
                'r_squared': 0.95
            }
        ]
        
        results = calculate_delta_s_prime(
            list_of_dicts,
            reference_cell_line_names="Reference",
            test_cell_line_names="Test"
        )
        
        assert "Reference" in results
        assert len(results["Reference"]) == 1
        assert results["Reference"][0]["delta_s_prime"] == -1.0  # 2.0 - 3.0
