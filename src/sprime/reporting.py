"""
Data quality reporting system for sprime.

Provides global configuration and reporting infrastructure for tracking
data quality issues, warnings, and processing metrics.
"""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Optional, List, Union
from datetime import datetime


class ConsoleOutput(Enum):
    """Console output verbosity levels."""
    NONE = "none"      # No console output
    SUMMARY = "summary"  # Brief summary (default)
    VERBOSE = "verbose"  # Detailed output with all warnings


class ReportingConfig:
    """Global configuration for data quality reporting."""
    
    # Log file settings
    log_to_file: bool = False
    log_filepath: Optional[Path] = None
    
    # Console output settings
    console_output: ConsoleOutput = ConsoleOutput.SUMMARY
    
    @classmethod
    def configure(
        cls,
        log_to_file: bool = False,
        log_filepath: Optional[Union[str, Path]] = None,
        console_output: Union[ConsoleOutput, str] = ConsoleOutput.SUMMARY
    ):
        """
        Configure global reporting settings.
        
        Args:
            log_to_file: If True, write detailed log file (default: False)
            log_filepath: Path to log file (default: auto-generated from input filename)
            console_output: Console verbosity - "none", "summary", or "verbose" (default: "summary")
        """
        cls.log_to_file = log_to_file
        
        if log_filepath:
            cls.log_filepath = Path(log_filepath)
        else:
            cls.log_filepath = None
        
        if isinstance(console_output, str):
            cls.console_output = ConsoleOutput(console_output.lower())
        else:
            cls.console_output = console_output
    
    @classmethod
    def reset(cls):
        """Reset to defaults."""
        cls.log_to_file = False
        cls.log_filepath = None
        cls.console_output = ConsoleOutput.SUMMARY


@dataclass
class WarningEntry:
    """Single warning with location information."""
    row_number: int  # CSV row number (1-based, header is row 1, 0 for non-CSV sources)
    category: str  # "DATA_QUALITY", "CURVE_FIT", "MISSING_DATA", "FORWARD_FILL", "NUMERICAL", "CALCULATION"
    message: str
    drug_id: Optional[str] = None
    compound_name: Optional[str] = None
    cell_line: Optional[str] = None
    field_name: Optional[str] = None
    
    def to_log_line(self) -> str:
        """Format for log file."""
        context_parts = []
        if self.drug_id:
            context_parts.append(f"Compound_ID: {self.drug_id}")
        if self.compound_name:
            context_parts.append(f"Compound: {self.compound_name}")
        if self.cell_line:
            context_parts.append(f"Cell_Line: {self.cell_line}")
        if self.field_name:
            context_parts.append(f"Field: {self.field_name}")
        
        context = " | ".join(context_parts) if context_parts else "N/A"
        row_str = f"Row {self.row_number:4d}" if self.row_number > 0 else "Row   N/A"
        return f"{row_str} | [{self.category:15s}] {self.message} | {context}"


@dataclass
class ProcessingReport:
    """Combined report for load + process operations."""
    
    # File info
    input_filepath: Optional[Path] = None
    
    # Summary metrics
    total_rows: int = 0
    rows_processed: int = 0
    rows_skipped: int = 0
    compounds_loaded: int = 0
    profiles_created: int = 0
    profiles_with_s_prime: int = 0
    profiles_failed_fit: int = 0
    
    # Data quality counts
    missing_drug_ids: int = 0
    missing_compound_names: int = 0
    missing_cell_lines: int = 0
    insufficient_data_points: int = 0
    invalid_numeric_values: int = 0
    forward_filled_fields: int = 0
    
    # All warnings (from both load and process)
    warnings: List[WarningEntry] = field(default_factory=list)
    
    def add_warning(self, row_number: int, category: str, message: str,
                   drug_id: Optional[str] = None,
                   compound_name: Optional[str] = None,
                   cell_line: Optional[str] = None,
                   field_name: Optional[str] = None):
        """Add a warning entry."""
        self.warnings.append(WarningEntry(
            row_number=row_number,
            category=category,
            message=message,
            drug_id=drug_id,
            compound_name=compound_name,
            cell_line=cell_line,
            field_name=field_name
        ))
    
    def write_log_file(self, filepath: Optional[Path] = None):
        """Write log file if enabled."""
        if not ReportingConfig.log_to_file:
            return
        
        # Use provided path or configured path or auto-generate
        if filepath:
            log_path = Path(filepath)
        elif ReportingConfig.log_filepath:
            log_path = ReportingConfig.log_filepath
        elif self.input_filepath:
            # Auto-generate from input filename
            log_path = self.input_filepath.parent / f"{self.input_filepath.stem}_processing.log"
        else:
            log_path = Path("sprime_processing.log")
        
        write_processing_log(self, log_path)
    
    def print_console_summary(self):
        """Print console summary based on configured verbosity."""
        if ReportingConfig.console_output == ConsoleOutput.NONE:
            return
        elif ReportingConfig.console_output == ConsoleOutput.SUMMARY:
            print_processing_summary(self)
        elif ReportingConfig.console_output == ConsoleOutput.VERBOSE:
            print_processing_summary_verbose(self)


def print_processing_summary(report: ProcessingReport):
    """Print brief summary to console."""
    print("\n" + "="*60)
    print("DATA PROCESSING SUMMARY")
    print("="*60)
    print(f"Total Rows:                 {report.total_rows}")
    print(f"Rows Processed:             {report.rows_processed}")
    print(f"Rows Skipped:                {report.rows_skipped}")
    print(f"Compounds Loaded:            {report.compounds_loaded}")
    print(f"Profiles Created:            {report.profiles_created}")
    print(f"Profiles with S' Calculated: {report.profiles_with_s_prime}")
    print(f"Profiles Failed:             {report.profiles_failed_fit}")
    print()
    
    if report.missing_drug_ids or report.missing_compound_names or \
       report.missing_cell_lines or report.insufficient_data_points or \
       report.invalid_numeric_values:
        print("DATA QUALITY ISSUES:")
        if report.missing_drug_ids:
            print(f"  Missing Compound IDs:       {report.missing_drug_ids}")
        if report.missing_compound_names:
            print(f"  Missing Compound Names:     {report.missing_compound_names}")
        if report.missing_cell_lines:
            print(f"  Missing Cell Lines:        {report.missing_cell_lines}")
        if report.insufficient_data_points:
            print(f"  Insufficient Data Points:   {report.insufficient_data_points}")
        if report.invalid_numeric_values:
            print(f"  Invalid Numeric Values:     {report.invalid_numeric_values}")
        print()
    
    if report.warnings:
        # Group by category
        by_category = {}
        for w in report.warnings:
            by_category.setdefault(w.category, []).append(w)
        
        print(f"WARNINGS: {len(report.warnings)} total")
        for category in sorted(by_category.keys()):
            warnings = by_category[category]
            print(f"  {category}: {len(warnings)}")
            
            # Show first 3 examples with row numbers
            for w in warnings[:3]:
                context = []
                if w.drug_id:
                    context.append(f"ID:{w.drug_id}")
                if w.cell_line:
                    context.append(f"CL:{w.cell_line}")
                context_str = f" ({', '.join(context)})" if context else ""
                row_str = f"Row {w.row_number}" if w.row_number > 0 else "N/A"
                print(f"    {row_str}: {w.message[:50]}{context_str}")
            
            if len(warnings) > 3:
                print(f"    ... and {len(warnings) - 3} more (see log file for details)")
        print()
    else:
        print("No warnings found.")
        print()
    
    print("="*60)


def print_processing_summary_verbose(report: ProcessingReport):
    """Print verbose summary with all warnings to console."""
    print("\n" + "="*60)
    print("DATA PROCESSING SUMMARY (VERBOSE)")
    print("="*60)
    print(f"Total Rows:                 {report.total_rows}")
    print(f"Rows Processed:             {report.rows_processed}")
    print(f"Rows Skipped:                {report.rows_skipped}")
    print(f"Compounds Loaded:            {report.compounds_loaded}")
    print(f"Profiles Created:            {report.profiles_created}")
    print(f"Profiles with S' Calculated: {report.profiles_with_s_prime}")
    print(f"Profiles Failed:             {report.profiles_failed_fit}")
    print()
    
    if report.missing_drug_ids or report.missing_compound_names or \
       report.missing_cell_lines or report.insufficient_data_points or \
       report.invalid_numeric_values:
        print("DATA QUALITY ISSUES:")
        if report.missing_drug_ids:
            print(f"  Missing Compound IDs:       {report.missing_drug_ids}")
        if report.missing_compound_names:
            print(f"  Missing Compound Names:     {report.missing_compound_names}")
        if report.missing_cell_lines:
            print(f"  Missing Cell Lines:        {report.missing_cell_lines}")
        if report.insufficient_data_points:
            print(f"  Insufficient Data Points:   {report.insufficient_data_points}")
        if report.invalid_numeric_values:
            print(f"  Invalid Numeric Values:     {report.invalid_numeric_values}")
        print()
    
    if report.warnings:
        # Group by category
        by_category = {}
        for w in report.warnings:
            by_category.setdefault(w.category, []).append(w)
        
        print(f"DETAILED WARNINGS: {len(report.warnings)} total\n")
        for category in sorted(by_category.keys()):
            warnings = by_category[category]
            print(f"[{category}] - {len(warnings)} warning(s)")
            print("-"*60)
            
            for warning in warnings:
                print(warning.to_log_line())
            print()
    else:
        print("No warnings found.")
        print()
    
    print("="*60)


def write_processing_log(report: ProcessingReport, log_filepath: Union[str, Path]):
    """
    Write detailed processing log to file.
    
    Args:
        report: ProcessingReport to write
        log_filepath: Path to output .log file
    """
    log_filepath = Path(log_filepath)
    
    with open(log_filepath, 'w', encoding='utf-8') as f:
        # Header
        f.write("="*80 + "\n")
        f.write("SPRIME DATA PROCESSING LOG\n")
        f.write("="*80 + "\n")
        if report.input_filepath:
            f.write(f"Input File: {report.input_filepath}\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n")
        
        # Summary
        f.write("SUMMARY\n")
        f.write("-"*80 + "\n")
        f.write(f"Total Rows in File:        {report.total_rows}\n")
        f.write(f"Rows Processed:            {report.rows_processed}\n")
        f.write(f"Rows Skipped:              {report.rows_skipped}\n")
        f.write(f"Compounds Loaded:          {report.compounds_loaded}\n")
        f.write(f"Profiles Created:          {report.profiles_created}\n")
        f.write(f"Profiles with S' Calculated: {report.profiles_with_s_prime}\n")
        f.write(f"Profiles Failed:           {report.profiles_failed_fit}\n")
        f.write("\n")
        
        # Data Quality Issues
        if report.missing_drug_ids or report.missing_compound_names or \
           report.missing_cell_lines or report.insufficient_data_points or \
           report.invalid_numeric_values:
            f.write("DATA QUALITY ISSUES\n")
            f.write("-"*80 + "\n")
            if report.missing_drug_ids:
                f.write(f"Missing Compound IDs:        {report.missing_drug_ids}\n")
            if report.missing_compound_names:
                f.write(f"Missing Compound Names:      {report.missing_compound_names}\n")
            if report.missing_cell_lines:
                f.write(f"Missing Cell Lines:         {report.missing_cell_lines}\n")
            if report.insufficient_data_points:
                f.write(f"Insufficient Data Points:   {report.insufficient_data_points}\n")
            if report.invalid_numeric_values:
                f.write(f"Invalid Numeric Values:     {report.invalid_numeric_values}\n")
            f.write("\n")
        
        # Warnings by category
        if report.warnings:
            f.write("DETAILED WARNINGS\n")
            f.write("-"*80 + "\n")
            
            # Group by category
            by_category = {}
            for w in report.warnings:
                by_category.setdefault(w.category, []).append(w)
            
            for category in sorted(by_category.keys()):
                warnings = by_category[category]
                f.write(f"\n[{category}] - {len(warnings)} warning(s)\n")
                f.write("-"*80 + "\n")
                
                for warning in warnings:
                    f.write(warning.to_log_line() + "\n")
        else:
            f.write("No warnings found.\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("End of log\n")
