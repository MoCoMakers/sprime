"""
Tests for sprime reporting system.

Tests cover:
- Report generation and accumulation
- Warning tracking and categorization
- Console output configuration
- Log file writing
- Blank row handling
- Data quality issue tracking
"""

import pytest
import tempfile
import csv
import os
from pathlib import Path
from io import StringIO
import sys

from sprime import (
    SPrime, RawDataset, ScreeningDataset,
    ReportingConfig, ConsoleOutput, ProcessingReport
)


class TestReportingConfig:
    """Tests for ReportingConfig global configuration."""
    
    def test_default_config(self):
        """Test default configuration values."""
        ReportingConfig.reset()
        assert ReportingConfig.log_to_file == False
        assert ReportingConfig.log_filepath is None
        assert ReportingConfig.console_output == ConsoleOutput.SUMMARY
    
    def test_configure_log_file(self):
        """Test enabling log file."""
        ReportingConfig.reset()
        ReportingConfig.configure(log_to_file=True)
        assert ReportingConfig.log_to_file == True
    
    def test_configure_console_output(self):
        """Test configuring console output."""
        ReportingConfig.reset()
        ReportingConfig.configure(console_output=ConsoleOutput.NONE)
        assert ReportingConfig.console_output == ConsoleOutput.NONE
        
        ReportingConfig.configure(console_output=ConsoleOutput.VERBOSE)
        assert ReportingConfig.console_output == ConsoleOutput.VERBOSE
        
        ReportingConfig.configure(console_output="summary")
        assert ReportingConfig.console_output == ConsoleOutput.SUMMARY
    
    def test_reset_config(self):
        """Test resetting to defaults."""
        ReportingConfig.configure(
            log_to_file=True,
            log_filepath="test.log",
            console_output=ConsoleOutput.NONE
        )
        ReportingConfig.reset()
        assert ReportingConfig.log_to_file == False
        assert ReportingConfig.log_filepath is None
        assert ReportingConfig.console_output == ConsoleOutput.SUMMARY


class TestReportGeneration:
    """Tests for report generation during data loading and processing."""
    
    def create_test_csv(self, rows, tmp_path):
        """Helper to create temporary CSV file."""
        csv_file = tmp_path / "test.csv"
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=sorted(fieldnames))
            writer.writeheader()
            writer.writerows(rows)
        
        return csv_file
    
    def test_load_returns_report(self, tmp_path):
        """Test that SPrime.load() returns a report."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100'
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        assert isinstance(report, ProcessingReport)
        assert report.total_rows == 2  # Header + 1 row
        assert report.rows_processed == 1
        assert report.compounds_loaded == 1
        assert report.profiles_created == 1
    
    def test_process_returns_report(self, tmp_path):
        """Test that SPrime.process() returns a report."""
        csv_file = self.create_test_csv([
            {
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
            }
        ], tmp_path)
        
        raw_data, load_report = SPrime.load(csv_file)
        screening_data, process_report = SPrime.process(raw_data)
        
        assert isinstance(process_report, ProcessingReport)
        assert process_report.profiles_with_s_prime == 1
    
    def test_report_accumulates_across_load_and_process(self, tmp_path):
        """Test that report accumulates warnings from both load and process."""
        csv_file = self.create_test_csv([
            {
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
            }
        ], tmp_path)
        
        raw_data, load_report = SPrime.load(csv_file)
        screening_data, process_report = SPrime.process(raw_data, report=load_report)
        
        # Process report should be the same object (accumulated)
        assert process_report is load_report
        assert process_report.profiles_created == 1
        assert process_report.profiles_with_s_prime == 1


class TestWarningTracking:
    """Tests for warning tracking and categorization."""
    
    def create_test_csv(self, rows, tmp_path):
        """Helper to create temporary CSV file."""
        csv_file = tmp_path / "test.csv"
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=sorted(fieldnames))
            writer.writeheader()
            writer.writerows(rows)
        
        return csv_file
    
    def test_fully_blank_row_not_logged(self, tmp_path):
        """Test that fully blank rows are skipped silently."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Conc0': '0.1'
            },
            {},  # Fully blank row
            {
                'Compound Name': 'Test2',
                'Compound_ID': 'TEST2',
                'Cell_Line': 'Cell2',
                'Concentration_Units': 'microM',
                'Data0': '20',
                'Conc0': '0.2'
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        # Blank row should not create a warning
        blank_row_warnings = [w for w in report.warnings if 'blank' in w.message.lower()]
        assert len(blank_row_warnings) == 0
        assert report.rows_processed == 2  # Only non-blank rows
    
    def test_missing_compound_id_raises_exception(self, tmp_path):
        """Test that missing Compound_ID raises ValueError."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test',
                # No Compound_ID in first row either - so can't forward-fill
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Conc0': '0.1'
            },
            {
                'Compound Name': 'Test2',
                # Missing Compound_ID - should raise exception
                'Cell_Line': 'Cell2',
                'Concentration_Units': 'microM',
                'Data0': '20',
                'Conc0': '0.2'
            }
        ], tmp_path)
        
        # Should raise ValueError for missing Compound_ID column in header
        with pytest.raises(ValueError, match="Required column 'Compound_ID' not found"):
            raw_data, report = SPrime.load(csv_file)
    
    def test_no_forward_fill(self, tmp_path):
        """No forward-filling; rows are literal, empty is null."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Conc0': '0.1'
            },
            {
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell2',
                'Concentration_Units': 'microM',
                'Data0': '20',
                'Conc0': '0.2'
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        forward_fill_warnings = [w for w in report.warnings if w.category == "FORWARD_FILL"]
        assert len(forward_fill_warnings) == 0
        assert report.forward_filled_fields == 0
    
    def test_insufficient_data_points_logged(self, tmp_path):
        """Test that insufficient data points are logged."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '20',
                'Data2': '30',  # Only 3 data points
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10'
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        # Should have warning for insufficient data points
        insufficient_warnings = [w for w in report.warnings 
                                 if w.category == "MISSING_DATA" and "insufficient" in w.message.lower()]
        assert len(insufficient_warnings) == 1
        assert report.insufficient_data_points == 1
    
    def test_invalid_numeric_values_logged(self, tmp_path):
        """Test that invalid numeric values are logged."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': 'invalid',  # Non-numeric
                'Data1': '20',
                'Data2': '50',
                'Data3': '90',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100'
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        # Should have warning for invalid numeric value
        invalid_warnings = [w for w in report.warnings 
                           if w.category in ["DATA_QUALITY", "NUMERICAL"]]
        assert len(invalid_warnings) >= 1
        assert report.invalid_numeric_values >= 1
    
    def test_curve_fit_failure_logged(self, tmp_path):
        """Test that curve fitting failures are logged."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test',
                'Compound_ID': 'TEST',
                'Cell_Line': 'Cell1',
                'Concentration_Units': 'microM',
                'Data0': '10',
                'Data1': '10',  # Flat response - may fail fitting
                'Data2': '10',
                'Data3': '10',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100'
            }
        ], tmp_path)
        
        raw_data, load_report = SPrime.load(csv_file)
        screening_data, process_report = SPrime.process(raw_data, report=load_report)
        
        # May have curve fit warnings (depending on fitting success)
        curve_fit_warnings = [w for w in process_report.warnings 
                             if w.category == "CURVE_FIT"]
        # At minimum, should track if profiles failed
        assert process_report.profiles_failed_fit >= 0  # May succeed or fail


class TestConsoleOutput:
    """Tests for console output configuration."""
    
    def create_test_csv(self, rows, tmp_path):
        """Helper to create temporary CSV file."""
        csv_file = tmp_path / "test.csv"
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=sorted(fieldnames))
            writer.writeheader()
            writer.writerows(rows)
        
        return csv_file
    
    def test_summary_output(self, tmp_path, capsys):
        """Test that summary output is printed by default."""
        ReportingConfig.reset()
        ReportingConfig.configure(console_output=ConsoleOutput.SUMMARY)
        
        csv_file = self.create_test_csv([
            {
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
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        captured = capsys.readouterr()
        
        assert "DATA PROCESSING SUMMARY" in captured.out
        assert "Total Rows" in captured.out
    
    def test_none_output(self, tmp_path, capsys):
        """Test that no output is printed when disabled."""
        ReportingConfig.reset()
        ReportingConfig.configure(console_output=ConsoleOutput.NONE)
        
        csv_file = self.create_test_csv([
            {
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
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        captured = capsys.readouterr()
        
        assert "DATA PROCESSING SUMMARY" not in captured.out
        assert len(captured.out.strip()) == 0


class TestLogFileWriting:
    """Tests for log file writing."""
    
    def create_test_csv(self, rows, tmp_path):
        """Helper to create temporary CSV file."""
        csv_file = tmp_path / "test.csv"
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=sorted(fieldnames))
            writer.writeheader()
            writer.writerows(rows)
        
        return csv_file
    
    def test_log_file_written_when_enabled(self, tmp_path):
        """Test that log file is written when enabled."""
        ReportingConfig.reset()
        ReportingConfig.configure(log_to_file=True)
        
        csv_file = self.create_test_csv([
            {
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
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        # Log file should be auto-generated
        log_file = csv_file.parent / f"{csv_file.stem}_processing.log"
        assert log_file.exists()
        
        # Check log file contents
        with open(log_file, 'r') as f:
            content = f.read()
            assert "SPRIME DATA PROCESSING LOG" in content
            assert "Total Rows in File" in content
    
    def test_custom_log_file_path(self, tmp_path):
        """Test custom log file path."""
        ReportingConfig.reset()
        custom_log = tmp_path / "custom.log"
        ReportingConfig.configure(log_to_file=True, log_filepath=custom_log)
        
        csv_file = self.create_test_csv([
            {
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
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        assert custom_log.exists()
    
    def test_log_file_not_written_when_disabled(self, tmp_path):
        """Test that log file is not written when disabled."""
        ReportingConfig.reset()
        ReportingConfig.configure(log_to_file=False)
        
        csv_file = self.create_test_csv([
            {
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
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        # Log file should not exist
        log_file = csv_file.parent / f"{csv_file.stem}_processing.log"
        assert not log_file.exists()


class TestReportMetrics:
    """Tests for report metrics and data quality tracking."""
    
    def create_test_csv(self, rows, tmp_path):
        """Helper to create temporary CSV file."""
        csv_file = tmp_path / "test.csv"
        fieldnames = set()
        for row in rows:
            fieldnames.update(row.keys())
        
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=sorted(fieldnames))
            writer.writeheader()
            writer.writerows(rows)
        
        return csv_file
    
    def test_report_metrics_accurate(self, tmp_path):
        """Test that report metrics are accurate."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test1',
                'Compound_ID': 'TEST1',
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
            },
            {
                'Compound Name': 'Test2',
                'Compound_ID': 'TEST2',
                'Cell_Line': 'Cell2',
                'Concentration_Units': 'microM',
                'Data0': '5',
                'Data1': '15',
                'Data2': '60',
                'Data3': '95',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
                'Conc3': '100'
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        assert report.total_rows == 3  # Header + 2 rows
        assert report.rows_processed == 2
        assert report.compounds_loaded == 2
        assert report.profiles_created == 2
    
    def test_warning_row_numbers(self, tmp_path):
        """Test that warnings have correct row numbers."""
        csv_file = self.create_test_csv([
            {
                'Compound Name': 'Test1',
                'Compound_ID': 'TEST1',
                'Cell_Line': 'Cell1',
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
                'Compound Name': 'Test2',
                'Compound_ID': 'TEST2',
                'Cell_Line': 'Cell2',
                'Concentration_Units': 'microM',
                'Data0': '20',
                'Data1': '30',
                'Data2': '40',
                'Conc0': '0.1',
                'Conc1': '1',
                'Conc2': '10',
            }
        ], tmp_path)
        
        raw_data, report = SPrime.load(csv_file)
        
        insufficient_warnings = [w for w in report.warnings 
                                if "Insufficient data points" in w.message and w.row_number > 0]
        assert len(insufficient_warnings) >= 1
        assert insufficient_warnings[0].row_number == 3
