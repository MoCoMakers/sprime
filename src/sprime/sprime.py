"""
sprime - A biomedical library for screening high-throughput screening data.

This module provides tools for analyzing dose-response curves and calculating
S' (S prime) values from quantitative high-throughput screening (QHTS) data.
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, List, Union, Tuple, Set, Literal
import csv
import math
import re
import warnings
from pathlib import Path

# Import hill fitting (may raise ImportError if scipy not installed)
try:
    from . import hill_fitting
except ImportError:
    hill_fitting = None  # Will be checked when fitting is attempted

# Import reporting (optional - for data quality tracking)
try:
    from .reporting import ProcessingReport, ReportingConfig
except ImportError:
    ProcessingReport = None
    ReportingConfig = None  # Fallback if reporting module not available


# ============================================================================
# Core Value Objects
# ============================================================================

@dataclass(frozen=True)
class Compound:
    """
    Represents a chemical compound/drug entity.
    
    Attributes:
        name: Compound name
        drug_id: Unique compound identifier (from Compound_ID column)
        pubchem_sid: PubChem substance identifier (optional)
        smiles: SMILES notation of chemical structure (optional)
    """
    name: str
    drug_id: str
    pubchem_sid: Optional[str] = None
    smiles: Optional[str] = None
    
    def __repr__(self):
        return f"Compound({self.name}, id={self.drug_id})"


@dataclass(frozen=True)
class CellLine:
    """
    Represents a cell line/clone.
    
    Attributes:
        name: Cell line name (e.g., 'ipNF96.11C', 'LS513_LARGE_INTESTINE')
        ref_id: Reference identifier (e.g., 'ACH-000007', 'depmap_id') (optional)
    """
    name: str
    ref_id: Optional[str] = None
    
    def __repr__(self):
        return f"CellLine({self.name})"


@dataclass
class Assay:
    """
    Represents a standardized biological assay procedure.
    
    Attributes:
        name: Assay name or identifier
        description: Assay description (optional)
        screen_id: Screen ID from data (e.g., 'HTS002') (optional)
        readout_type: Type of measurement (e.g., 'activity', 'luminescence') (optional)
        time_profile: Time point if applicable (e.g., '24Hr', '48Hr', '4Day') (optional)
    """
    name: str
    description: Optional[str] = None
    screen_id: Optional[str] = None
    readout_type: Optional[str] = None
    time_profile: Optional[str] = None
    
    def __repr__(self):
        return f"Assay({self.name})"


# ============================================================================
# Hill Curve Parameters
# ============================================================================

@dataclass
class HillCurveParams:
    """
    Parameters from four-parameter Hill equation fit.
    
    Hill equation: y = D + (A - D) / (1 + (x/C)^n)
    Where:
        A = upper (upper asymptote)
        D = lower (lower asymptote)
        C = ec50 (half-maximal concentration)
        n = hill_coefficient (slope/steepness)
    
    Attributes:
        ec50: Half-maximal effect concentration (AC50/EC50)
        upper: Upper asymptote (A parameter)
        lower: Lower asymptote (D parameter)
        hill_coefficient: Hill coefficient/slope (n parameter) (optional)
        r_squared: R-squared goodness of fit statistic (optional)
    """
    ec50: float
    upper: float
    lower: float
    hill_coefficient: Optional[float] = None
    r_squared: Optional[float] = None
    
    @property
    def amplitude(self) -> float:
        """Response amplitude (Upper - Lower)"""
        return self.upper - self.lower


# ============================================================================
# Dose-Response Profile
# ============================================================================

@dataclass
class DoseResponseProfile:
    """
    Represents a dose-response profile for one Compound-CellLine pair in one Assay.
    
    Contains raw data, fitted curve parameters, and calculated S' value.
    
    Attributes:
        compound: Compound entity
        cell_line: CellLine entity
        assay: Assay entity
        concentrations: Raw concentration values (list of floats) (optional)
        responses: Raw response values (list of floats) (optional)
        concentration_units: Units for concentrations (default: 'microM')
        hill_params: Fitted Hill curve parameters (optional)
        s_prime: Calculated S' value (optional)
        rank: Rank of S' value (optional)
        metadata: Additional metadata from CSV (e.g., MOA, drug targets) (optional)
    """
    compound: Compound
    cell_line: CellLine
    assay: Assay
    
    # Raw data
    concentrations: Optional[List[float]] = None
    responses: Optional[List[float]] = None
    concentration_units: str = "microM"
    
    # Fitted parameters
    hill_params: Optional[HillCurveParams] = None
    
    # Results
    s_prime: Optional[float] = None
    rank: Optional[int] = None
    
    # Additional metadata
    metadata: Optional[Dict[str, str]] = None
    
    def fit_hill_curve(self, **fit_params) -> HillCurveParams:
        """
        Fit four-parameter Hill equation to raw dose-response data.
        
        Updates self.hill_params with fitted parameters.
        
        Args:
            **fit_params: Additional parameters for curve fitting
            
        Returns:
            HillCurveParams: Fitted curve parameters
            
        Raises:
            ValueError: If raw data is not available or hill_params already exist
        """
        if self.hill_params is not None:
            return self.hill_params  # Already fitted
        
        if self.concentrations is None or self.responses is None:
            raise ValueError("Need raw data to fit Hill curve")
        
        if len(self.concentrations) != len(self.responses):
            raise ValueError("Concentrations and responses must have same length")
        
        if hill_fitting is None:
            raise ImportError(
                "Hill curve fitting requires scipy. "
                "Install with: pip install scipy"
            )
        
        # Fit Hill curve using hill_fitting module
        self.hill_params = hill_fitting.fit_hill_curve(
            self.concentrations,
            self.responses,
            **fit_params
        )
        
        return self.hill_params
    
    def calculate_s_prime(self) -> float:
        """
        Calculate S' = asinh((Upper-Lower)/EC50).
        
        Requires hill_params to be set (call fit_hill_curve first).
        
        Returns:
            float: S' value
            
        Raises:
            ValueError: If hill_params is not available
        """
        if self.hill_params is None:
            raise ValueError("Must fit Hill curve before calculating S'")
        
        # S' = asinh((A-D)/C) where A=upper, D=lower, C=EC50
        ratio = self.hill_params.amplitude / self.hill_params.ec50
        self.s_prime = math.asinh(ratio)  # asinh(x) = ln(x + sqrt(x^2 + 1))
        return self.s_prime
    
    def fit_and_calculate_s_prime(self, **fit_params) -> float:
        """
        Convenience method: fit curve then calculate S'.
        
        Args:
            **fit_params: Parameters for curve fitting
            
        Returns:
            float: S' value
        """
        self.fit_hill_curve(**fit_params)
        return self.calculate_s_prime()


# ============================================================================
# Validation Helpers
# ============================================================================

def _validate_required_columns(
    columns: List[str],
    source_name: str = "data",
    values_as: str = "columns",
) -> None:
    """
    Validate that required columns exist in the data structure.
    
    Raises ValueError with clear error messages if required columns are missing.
    
    Args:
        columns: List of column names (from CSV header or dict keys)
        source_name: Name of the data source (for error messages)
        values_as: "columns" (DATA*/CONC*) or "list" (Responses, Concentrations)
        
    Raises:
        ValueError: If required columns are missing
    """
    columns_lower = {col.lower(): col for col in columns}  # Case-insensitive lookup
    columns_set = set(columns)
    
    # Check for Cell_Line (required)
    if 'Cell_Line' not in columns_set and 'cell_line' not in columns_lower:
        available = ', '.join(sorted(columns)) if columns else '(none)'
        raise ValueError(
            f"Required column 'Cell_Line' not found in {source_name}. "
            f"Available columns: {available}"
        )
    
    # Check for Compound_ID (required). NCGCID is optional pass-through.
    has_compound_id = 'Compound_ID' in columns_set or 'compound_id' in columns_lower
    
    if not has_compound_id:
        available = ', '.join(sorted(columns)) if columns else '(none)'
        raise ValueError(
            f"Required column 'Compound_ID' not found in {source_name}. "
            f"Available columns: {available}"
        )
    
    # Raw data: columns (DATA*/CONC*) or list (Responses, Concentrations)
    has_raw_columns = False
    has_raw_list = False
    if values_as == "list":
        has_raw_list = (
            ('Responses' in columns_set or 'responses' in columns_lower) and
            ('Concentrations' in columns_set or 'concentrations' in columns_lower)
        )
    else:
        has_raw_data = any(
            col.startswith('Data') or col.startswith('DATA') or
            col.startswith('data') or col.startswith('Data')
            for col in columns
        )
        has_raw_conc = any(
            (col.startswith('Conc') or col.startswith('CONC') or
             col.startswith('conc') or col.startswith('Conc'))
            and 'Units' not in col and 'units' not in col
            for col in columns
        )
        has_raw_columns = has_raw_data and has_raw_conc
    
    has_raw_data_pattern = has_raw_columns or has_raw_list
    has_ac50 = 'AC50' in columns_set or 'ac50' in columns_lower or 'ec50' in columns_lower
    has_upper = 'Upper' in columns_set or 'upper' in columns_lower or 'Infinity' in columns_set or 'infinity' in columns_lower
    has_lower = 'Lower' in columns_set or 'lower' in columns_lower or 'Zero' in columns_set or 'zero' in columns_lower
    has_precalc_params = has_ac50 and has_upper and has_lower
    
    if not (has_raw_data_pattern or has_precalc_params):
        available = ', '.join(sorted(columns)) if columns else '(none)'
        raise ValueError(
            f"No dose-response data found in {source_name}. "
            f"Expected either: (1) raw data columns (DATA*/CONC* or Responses/Concentrations), or "
            f"(2) pre-calculated parameters (AC50/Upper/Lower or ec50/Infinity/Zero). "
            f"Available columns: {available}"
        )
    
    # Concentration_Units required when raw data present (Path A)
    has_concentration_units = 'Concentration_Units' in columns_set or 'concentration_units' in columns_lower
    if has_raw_data_pattern and not has_concentration_units:
        available = ', '.join(sorted(columns)) if columns else '(none)'
        raise ValueError(
            f"Required column 'Concentration_Units' not found in {source_name}. "
            f"Raw dose-response data requires concentration units. "
            f"Available columns: {available}"
        )


def _reserved_column_names(values_as: str, columns: List[str]) -> Set[str]:
    """
    Return the set of reserved column names (required + structural + pre-calc only).
    Used to exclude these from generic metadata pass-through at S' load. All other
    columns are stored as metadata under exact header, value as-is.
    """
    reserved = {
        "Compound_ID", "Cell_Line", "Compound Name", "pubchem_sid", "SMILES", "NCGCID",
        "Cell_Line_Ref_ID", "Concentration_Units",
        "AC50", "ec50", "Upper", "Lower", "Infinity", "Zero", "Hill_Slope", "Hill", "slope",
        "r2", "R²", "Rank", "S'", "S Prime", "Responses", "Concentrations",
    }
    if values_as == "columns":
        for c in columns:
            if re.match(r"^data\d+$", c, re.IGNORECASE) or re.match(r"^conc\d+$", c, re.IGNORECASE):
                reserved.add(c)
    return reserved


def _resolve_moa(meta: Optional[Dict]) -> str:
    """Resolve MOA from metadata (checks MOA, MoA, moa). Used only for delta S' output."""
    if not meta:
        return ""
    for k in ("MOA", "MoA", "moa"):
        v = meta.get(k, "")
        v = v if isinstance(v, str) else str(v)
        if v.strip():
            return v
    return ""


def _resolve_drug_targets(meta: Optional[Dict]) -> str:
    """Resolve drug targets from metadata (checks drug targets, Target, target). Used only for delta S' output."""
    if not meta:
        return ""
    for k in ("drug targets", "Target", "target"):
        v = meta.get(k, "")
        v = v if isinstance(v, str) else str(v)
        if v.strip():
            return v
    return ""


def _convert_dataframe_to_dict_list(df) -> List[Dict]:
    """
    Convert pandas DataFrame to list of dictionaries.
    
    Args:
        df: pandas DataFrame
        
    Returns:
        List of dictionaries (one per row)
        
    Raises:
        ImportError: If pandas is not installed
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError(
            "pandas is required to load DataFrames. Install with: pip install pandas"
        )
    
    if not isinstance(df, pd.DataFrame):
        raise TypeError(f"Expected pandas DataFrame, got {type(df).__name__}")
    
    return df.to_dict('records')


# ============================================================================
# Raw Dataset - Loads and holds raw data
# ============================================================================

class RawDataset:
    """
    Holds raw dose-response curve data loaded from CSV/files.
    
    Focus: Data loading, validation, basic structure.
    May contain profiles with raw data (concentrations/responses) or 
    pre-calculated Hill params from CSV.
    """
    
    def __init__(self, assay: Assay):
        """
        Initialize RawDataset.
        
        Args:
            assay: Assay entity for this dataset
        """
        self.assay = assay
        self._profiles: Dict[Tuple[str, str], DoseResponseProfile] = {}  # (compound_id, cellline_name) -> profile
    
    @classmethod
    def load_from_file(
        cls,
        filepath: Union[str, Path],
        assay_name: Optional[str] = None,
        report: Optional['ProcessingReport'] = None,
        values_as: str = "columns",
        **csv_kwargs
    ) -> Tuple['RawDataset', 'ProcessingReport']:
        """
        Load CSV file and create RawDataset with quality reporting.
        
        Handles:
        - Extracting raw data (Data0..DataN, Conc0..ConcN) or Responses/Concentrations
        - Loading pre-calculated params (AC50, Upper, Lower) if present
        - Creating DoseResponseProfile objects
        - Tracking data quality issues and warnings
        
        Args:
            filepath: Path to CSV file
            assay_name: Name for assay (defaults to filename stem)
            report: Optional ProcessingReport to accumulate warnings (creates new if None)
            values_as: "columns" (DATA*/CONC*) or "list" (Responses, Concentrations)
            **csv_kwargs: Additional arguments for csv.DictReader
            
        Returns:
            Tuple of (RawDataset, ProcessingReport)
        """
        filepath = Path(filepath)
        
        # Create or use existing report
        if report is None:
            if ProcessingReport is not None:
                report = ProcessingReport()
                report.input_filepath = filepath
            else:
                # Fallback if reporting not available
                report = None
        
        # Infer assay name
        if assay_name is None:
            assay_name = filepath.stem
        
        assay = Assay(name=assay_name)
        raw_dataset = cls(assay=assay)
        
        # Read CSV
        with open(filepath, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, **csv_kwargs)
            fieldnames = reader.fieldnames
            rows = list(reader)
        
        # Validate required columns exist in header
        if fieldnames:
            _validate_required_columns(
                fieldnames, source_name=f"CSV file '{filepath}'", values_as=values_as
            )
        
        if report:
            report.total_rows = len(rows) + 1  # +1 for header
        
        if not rows:
            if report is None:
                # Return dummy report if reporting not available
                return raw_dataset, None
            return raw_dataset, report
        
        reserved = _reserved_column_names(values_as, fieldnames or [])
        
        # Track compounds seen
        compounds_seen = set()
        
        # Process rows with line numbers (row 1 is header, so start at 2)
        for row_num, row in enumerate(rows, start=2):
            # Check if row is fully blank (all values empty/whitespace)
            is_fully_blank = not any(
                v.strip() if isinstance(v, str) else str(v).strip()
                for v in row.values() if v
            )
            
            if is_fully_blank:
                # Skip fully blank rows silently (no logging, don't count as processed)
                continue
            
            if report:
                report.rows_processed += 1
            
            # Check for empty cell line - RAISE EXCEPTION
            cell_line_name = row.get('Cell_Line', '').strip()
            if not cell_line_name:
                raise ValueError(
                    f"Row {row_num}: Missing required 'Cell_Line' value. "
                    f"All rows must have a cell line specified."
                )
            
            # Get compound info (Compound_ID required; NCGCID is pass-through only)
            # Rows are taken literally: empty values are null, no forward-filling.
            compound_name = row.get('Compound Name', '').strip() or 'Unknown'
            compound_id = row.get('Compound_ID', '').strip()
            
            if not compound_id:
                raise ValueError(
                    f"Row {row_num}: Missing required 'Compound_ID' value. "
                    f"Compound: {compound_name}, Cell_Line: {cell_line_name}. "
                    f"All rows must have a compound identifier."
                )
            
            # Track compound
            if report and compound_id not in compounds_seen:
                compounds_seen.add(compound_id)
                report.compounds_loaded += 1
            
            # Check for missing compound name
            if report and compound_name == 'Unknown':
                report.add_warning(
                    row_number=row_num,
                    category="DATA_QUALITY",
                    message="Compound Name missing, using 'Unknown'",
                    drug_id=compound_id,
                    cell_line=cell_line_name,
                    field_name="Compound Name"
                )
                report.missing_compound_names += 1
            
            compound = Compound(
                name=compound_name,
                drug_id=compound_id,
                pubchem_sid=row.get('pubchem_sid', '').strip() or None,
                smiles=row.get('SMILES', '').strip() or None
            )
            
            # Create CellLine
            cell_line = CellLine(
                name=cell_line_name,
                ref_id=row.get('Cell_Line_Ref_ID', '').strip() or None
            )
            
            # Extract raw dose-response data (if present)
            concentrations = None
            responses = None
            
            if values_as == "list":
                resp_key = next((k for k in row if k.lower() == 'responses'), None)
                conc_key = next((k for k in row if k.lower() == 'concentrations'), None)
                if resp_key and conc_key:
                    resp_str = (row.get(resp_key) or '').strip()
                    conc_str = (row.get(conc_key) or '').strip()
                    if resp_str and conc_str:
                        responses = []
                        concentrations = []
                        for part in resp_str.split(','):
                            t = part.strip()
                            if t:
                                try:
                                    v = float(t)
                                    if report and (math.isnan(v) or math.isinf(v)):
                                        report.add_warning(
                                            row_number=row_num,
                                            category="NUMERICAL",
                                            message=f"Invalid numeric value (NaN/Inf) in Responses",
                                            drug_id=compound_id,
                                            compound_name=compound_name,
                                            cell_line=cell_line_name,
                                            field_name="Responses",
                                        )
                                        report.invalid_numeric_values += 1
                                        continue
                                    responses.append(v)
                                except (ValueError, TypeError) as e:
                                    if report:
                                        report.add_warning(
                                            row_number=row_num,
                                            category="DATA_QUALITY",
                                            message=f"Non-numeric value in Responses: {e}",
                                            drug_id=compound_id,
                                            compound_name=compound_name,
                                            cell_line=cell_line_name,
                                            field_name="Responses",
                                        )
                                        report.invalid_numeric_values += 1
                                    continue
                        conc_parts = []
                        for part in conc_str.split(','):
                            t = part.strip()
                            if t:
                                try:
                                    v = float(t)
                                    if report and (math.isnan(v) or math.isinf(v)):
                                        report.add_warning(
                                            row_number=row_num,
                                            category="NUMERICAL",
                                            message=f"Invalid numeric value (NaN/Inf) in Concentrations",
                                            drug_id=compound_id,
                                            compound_name=compound_name,
                                            cell_line=cell_line_name,
                                            field_name="Concentrations",
                                        )
                                        report.invalid_numeric_values += 1
                                        continue
                                    conc_parts.append(v)
                                except (ValueError, TypeError) as e:
                                    if report:
                                        report.add_warning(
                                            row_number=row_num,
                                            category="DATA_QUALITY",
                                            message=f"Non-numeric value in Concentrations: {e}",
                                            drug_id=compound_id,
                                            compound_name=compound_name,
                                            cell_line=cell_line_name,
                                            field_name="Concentrations",
                                        )
                                        report.invalid_numeric_values += 1
                                    continue
                        if len(responses) != len(conc_parts):
                            raise ValueError(
                                f"Row {row_num}: Responses and Concentrations length mismatch "
                                f"({len(responses)} vs {len(conc_parts)}). "
                                f"Compound: {compound_name}, Cell_Line: {cell_line_name}."
                            )
                        if len(responses) < 4:
                            if report:
                                report.add_warning(
                                    row_number=row_num,
                                    category="MISSING_DATA",
                                    message=f"Insufficient data points: {len(responses)} found, need 4+ for curve fitting",
                                    drug_id=compound_id,
                                    compound_name=compound_name,
                                    cell_line=cell_line_name,
                                )
                                report.insufficient_data_points += 1
                            concentrations = None
                            responses = None
                        else:
                            units = row.get('Concentration_Units', '').strip()
                            if not units:
                                raise ValueError(
                                    f"Row {row_num}: Missing required 'Concentration_Units' for raw data. "
                                    f"Compound: {compound_name}, Cell_Line: {cell_line_name}. "
                                    f"Raw dose-response data requires Concentration_Units."
                                )
                            concentrations = convert_to_micromolar(conc_parts, units)
                            responses = responses
            
            elif True:
                # Columns format: DATA*/CONC*
                data_cols = [k for k in row.keys() if k.startswith('Data') or k.startswith('DATA')]
                conc_cols = [k for k in row.keys() 
                            if (k.startswith('Conc') or k.startswith('CONC')) 
                            and 'Units' not in k and 'units' not in k]
                if data_cols and conc_cols:
                    # Sort columns to ensure correct order
                    data_cols = sorted(data_cols, key=lambda x: int(''.join(filter(str.isdigit, x)) or '0'))
                    conc_cols = sorted(conc_cols, key=lambda x: int(''.join(filter(str.isdigit, x)) or '0'))
                    
                    responses = []
                    concentrations = []
                    
                    for data_col, conc_col in zip(data_cols, conc_cols):
                        try:
                            resp_val = row.get(data_col, '') or ''
                            conc_val = row.get(conc_col, '') or ''
                            resp_val = resp_val.strip() if isinstance(resp_val, str) else ''
                            conc_val = conc_val.strip() if isinstance(conc_val, str) else ''
                            if resp_val and conc_val:
                                resp_float = float(resp_val)
                                conc_float = float(conc_val)
                                
                                # Check for invalid numeric values
                                if report:
                                    if math.isnan(resp_float) or math.isinf(resp_float):
                                        report.add_warning(
                                            row_number=row_num,
                                            category="NUMERICAL",
                                            message=f"Invalid numeric value (NaN/Inf) in {data_col}",
                                            drug_id=compound_id,
                                            compound_name=compound_name,
                                            cell_line=cell_line_name,
                                            field_name=data_col
                                        )
                                        report.invalid_numeric_values += 1
                                        continue
                                    
                                    if math.isnan(conc_float) or math.isinf(conc_float):
                                        report.add_warning(
                                            row_number=row_num,
                                            category="NUMERICAL",
                                            message=f"Invalid numeric value (NaN/Inf) in {conc_col}",
                                            drug_id=compound_id,
                                            compound_name=compound_name,
                                            cell_line=cell_line_name,
                                            field_name=conc_col
                                        )
                                        report.invalid_numeric_values += 1
                                        continue
                                
                                responses.append(resp_float)
                                concentrations.append(conc_float)
                        except (ValueError, TypeError) as e:
                            if report:
                                report.add_warning(
                                    row_number=row_num,
                                    category="DATA_QUALITY",
                                    message=f"Non-numeric value in {data_col} or {conc_col}: {str(e)}",
                                    drug_id=compound_id,
                                    compound_name=compound_name,
                                    cell_line=cell_line_name,
                                    field_name=f"{data_col}/{conc_col}"
                                )
                                report.invalid_numeric_values += 1
                            continue
                    
                    # Check for sufficient data points
                    if report and len(responses) < 4:
                        report.add_warning(
                            row_number=row_num,
                            category="MISSING_DATA",
                            message=f"Insufficient data points: {len(responses)} found, need 4+ for curve fitting",
                            drug_id=compound_id,
                            compound_name=compound_name,
                            cell_line=cell_line_name
                        )
                        report.insufficient_data_points += 1
                    
                    if concentrations and responses:
                        units = row.get('Concentration_Units', '').strip()
                        if not units:
                            raise ValueError(
                                f"Row {row_num}: Missing required 'Concentration_Units' for raw data. "
                                f"Compound: {compound_name}, Cell_Line: {cell_line_name}. "
                                f"Raw dose-response data requires Concentration_Units."
                            )
                        concentrations = convert_to_micromolar(concentrations, units)
                        concentrations = concentrations if concentrations else None
                        responses = responses if responses else None
            
            # Extract pre-calculated Hill params (if present)
            hill_params = None
            ac50 = row.get('AC50', '').strip() or row.get('ec50', '').strip()
            if ac50:
                try:
                    r2_raw = row.get('r2', row.get('R²', ''))
                    r_sq = _try_float(r2_raw)
                    hill_params = HillCurveParams(
                        ec50=float(ac50),
                        upper=float(row.get('Upper', row.get('Infinity', '0')).strip() or '0'),
                        lower=float(row.get('Lower', row.get('Zero', '0')).strip() or '0'),
                        hill_coefficient=_try_float(row.get('Hill_Slope', row.get('Hill', row.get('slope', '')))),
                        r_squared=r_sq
                    )
                except (ValueError, TypeError):
                    hill_params = None
            
            # Pre-calculated S' (if present)
            s_prime = _try_float(row.get("S'", row.get('S Prime', '')))
            rank = _try_int(row.get('Rank', ''))
            
            # Validate that row has either raw data or pre-calculated params
            has_raw_data = concentrations is not None and responses is not None and len(concentrations) > 0
            has_precalc_params = hill_params is not None

            if not (has_raw_data or has_precalc_params):
                raise ValueError(
                    f"Row {row_num}: No dose-response data found for compound '{compound_name}' "
                    f"(Compound_ID: {compound_id}) in cell line '{cell_line_name}'. "
                    f"Row must have either: (1) raw data columns (DATA*/CONC*), or "
                    f"(2) pre-calculated parameters (AC50/Upper/Lower)."
                )
            
            # Extract metadata: generic pass-through (all non-reserved columns, exact header, value as-is)
            metadata = {}
            for col in row.keys():
                if col in reserved:
                    continue
                raw = row.get(col, "")
                raw = raw if isinstance(raw, str) else str(raw)
                metadata[col] = raw
            
            # Create profile
            profile = DoseResponseProfile(
                compound=compound,
                cell_line=cell_line,
                assay=assay,
                concentrations=concentrations,
                responses=responses,
                concentration_units="microM",
                hill_params=hill_params,
                s_prime=s_prime,
                rank=rank,
                metadata=metadata if metadata else None
            )
            
            raw_dataset.add_profile(profile)
            if report:
                report.profiles_created += 1
        
        if report is None:
            # Return dummy report if reporting not available
            return raw_dataset, None
        return raw_dataset, report
    
    def add_profile(self, profile: DoseResponseProfile):
        """
        Add a profile to the dataset.
        
        Args:
            profile: DoseResponseProfile to add
            
        Raises:
            ValueError: If profile for this compound-cellline pair already exists
        """
        key = (profile.compound.drug_id, profile.cell_line.name)
        if key in self._profiles:
            raise ValueError(f"Profile for {key} already exists")
        self._profiles[key] = profile
    
    def get_profile(self, compound: Union[Compound, str], 
                    cell_line: Union[CellLine, str]) -> Optional[DoseResponseProfile]:
        """
        Retrieve a specific profile.
        
        Args:
            compound: Compound object or drug_id string
            cell_line: CellLine object or cell_line name string
            
        Returns:
            DoseResponseProfile or None if not found
        """
        compound_id = compound.drug_id if isinstance(compound, Compound) else compound
        cellline_name = cell_line.name if isinstance(cell_line, CellLine) else cell_line
        return self._profiles.get((compound_id, cellline_name))
    
    def to_screening_dataset(
        self,
        report: Optional['ProcessingReport'] = None,
        allow_overwrite_hill_coefficients: bool = False,
        **fit_params
    ) -> Tuple['ScreeningDataset', 'ProcessingReport']:
        """
        Process raw data into ScreeningDataset with quality reporting.
        
        For each profile:
        1. Fit Hill curve when raw data exists, or use pre-calc when no raw.
        2. Always compute S' from Hill params; warn if S' was provided in input and is overwritten.
        
        When both raw data and pre-calc (AC50/Upper/Lower/Hill_Slope/r2) exist:
        - allow_overwrite_hill_coefficients=False (default): raise (would overwrite).
        - allow_overwrite_hill_coefficients=True: fit, overwrite pre-calc, and log
          a warning that pre-calc Hill params were overwritten.
        
        Args:
            report: Optional ProcessingReport to accumulate warnings (creates new if None)
            allow_overwrite_hill_coefficients: If True, allow overwriting pre-calc
                Hill params (AC50, Upper, Lower, Hill_Slope, r2) with fitted values.
                Default False: raise when we would overwrite.
            **fit_params: Parameters for curve fitting (e.g. maxfev, bounds).
            
        Returns:
            Tuple of (ScreeningDataset, ProcessingReport)
            
        Raises:
            ValueError: If profile cannot be processed, or would overwrite
                pre-calc Hill coefficients without allow_overwrite_hill_coefficients.
        """
        if report is None:
            if ProcessingReport is not None:
                report = ProcessingReport()
            else:
                report = None
        
        screening_dataset = ScreeningDataset(assay=self.assay)
        
        for profile in self._profiles.values():
            processed_profile = DoseResponseProfile(
                compound=profile.compound,
                cell_line=profile.cell_line,
                assay=profile.assay,
                concentrations=profile.concentrations,
                responses=profile.responses,
                concentration_units=profile.concentration_units,
                hill_params=profile.hill_params,
                s_prime=profile.s_prime,
                rank=profile.rank,
                metadata=profile.metadata
            )
            
            has_raw = (
                processed_profile.concentrations is not None
                and processed_profile.responses is not None
                and len(processed_profile.concentrations) > 0
            )
            had_precalc = processed_profile.hill_params is not None
            
            if has_raw:
                if had_precalc and not allow_overwrite_hill_coefficients:
                    raise ValueError(
                        f"Pre-calculated Hill parameters (AC50, Upper, Lower, Hill_Slope, r2) would be "
                        f"overwritten by fitted values for compound '{processed_profile.compound.name}' "
                        f"(Compound_ID: {processed_profile.compound.drug_id}) in cell line "
                        f"'{processed_profile.cell_line.name}'. "
                        f"Set allow_overwrite_hill_coefficients=True to permit."
                    )
                processed_profile.hill_params = None
                try:
                    processed_profile.fit_hill_curve(**fit_params)
                except (RuntimeError, ValueError) as e:
                    raise ValueError(
                        f"Curve fitting failed for compound '{processed_profile.compound.name}' "
                        f"(Compound_ID: {processed_profile.compound.drug_id}) in cell line "
                        f"'{processed_profile.cell_line.name}': {str(e)}"
                    )
                if had_precalc:
                    msg = (
                        f"Pre-calc Hill parameters (AC50, Upper, Lower, Hill_Slope, r2) overwritten by "
                        f"fitted values for '{processed_profile.compound.name}' / "
                        f"'{processed_profile.cell_line.name}'."
                    )
                    warnings.warn(msg, UserWarning, stacklevel=2)
                    if report:
                        report.add_warning(
                            row_number=0,
                            category="OVERWRITE_HILL",
                            message=msg,
                            drug_id=processed_profile.compound.drug_id,
                            compound_name=processed_profile.compound.name,
                            cell_line=processed_profile.cell_line.name,
                        )
            else:
                if processed_profile.hill_params is None:
                    raise ValueError(
                        f"No data available to process for compound '{processed_profile.compound.name}' "
                        f"(Compound_ID: {processed_profile.compound.drug_id}) in cell line "
                        f"'{processed_profile.cell_line.name}'. Profile has neither raw data "
                        f"(concentrations/responses) nor pre-calculated Hill curve parameters."
                    )
                msg = (
                    f"Using pre-calc Hill parameters as-is for '{processed_profile.compound.name}' / "
                    f"'{processed_profile.cell_line.name}' (no raw data)."
                )
                warnings.warn(msg, UserWarning, stacklevel=2)
                if report:
                    report.add_warning(
                        row_number=0,
                        category="PRE_CALC_ONLY",
                        message=msg,
                        drug_id=processed_profile.compound.drug_id,
                        compound_name=processed_profile.compound.name,
                        cell_line=processed_profile.cell_line.name,
                    )
            
            # Check fit quality
            if report and processed_profile.hill_params and processed_profile.hill_params.r_squared is not None:
                if processed_profile.hill_params.r_squared < 0.7:
                    report.add_warning(
                        row_number=0,
                        category="CURVE_FIT",
                        message=f"Poor fit quality: R² = {processed_profile.hill_params.r_squared:.3f}",
                        drug_id=processed_profile.compound.drug_id,
                        compound_name=processed_profile.compound.name,
                        cell_line=processed_profile.cell_line.name
                    )
            
            # Always compute S'; warn if overwriting existing S' from CSV
            if processed_profile.s_prime is not None:
                msg = (
                    f"S' overwritten by recomputation for '{processed_profile.compound.name}' / "
                    f"'{processed_profile.cell_line.name}' (S' was provided in input)."
                )
                warnings.warn(msg, UserWarning, stacklevel=2)
                if report:
                    report.add_warning(
                        row_number=0,
                        category="OVERWRITE_S_PRIME",
                        message=msg,
                        drug_id=processed_profile.compound.drug_id,
                        compound_name=processed_profile.compound.name,
                        cell_line=processed_profile.cell_line.name,
                    )
            try:
                processed_profile.calculate_s_prime()
            except Exception as e:
                raise ValueError(
                    f"S' calculation failed for compound '{processed_profile.compound.name}' "
                    f"(Compound_ID: {processed_profile.compound.drug_id}) in cell line "
                    f"'{processed_profile.cell_line.name}': {str(e)}"
                )
            
            # Validate S' was calculated
            if processed_profile.s_prime is None:
                raise ValueError(
                    f"S' calculation returned None for compound '{processed_profile.compound.name}' "
                    f"(Compound_ID: {processed_profile.compound.drug_id}) in cell line "
                    f"'{processed_profile.cell_line.name}'. This indicates invalid Hill curve parameters."
                )
            
            # Add profile to screening dataset
            try:
                screening_dataset.add_profile(processed_profile)
                if report:
                    report.profiles_with_s_prime += 1
            except ValueError as e:
                # Profile already exists or invalid - re-raise with context
                raise ValueError(
                    f"Cannot add profile for compound '{processed_profile.compound.name}' "
                    f"(Compound_ID: {processed_profile.compound.drug_id}) in cell line "
                    f"'{processed_profile.cell_line.name}': {str(e)}"
                )
        
        if report is None:
            # Return dummy report if reporting not available
            return screening_dataset, None
        return screening_dataset, report
    
    @property
    def profiles(self):
        """Iterator over all profiles"""
        return self._profiles.values()
    
    def __len__(self):
        return len(self._profiles)


# ============================================================================
# Screening Dataset - Processed data ready for analysis
# ============================================================================

class ScreeningDataset:
    """
    Processed dataset with all S' values calculated.
    
    All profiles in this dataset:
    - Have fitted Hill curve parameters
    - Have calculated S' values
    - Are ready for analysis operations
    """
    
    def __init__(self, assay: Assay):
        """
        Initialize ScreeningDataset.
        
        Args:
            assay: Assay entity for this dataset
        """
        self.assay = assay
        self._profiles: Dict[Tuple[str, str], DoseResponseProfile] = {}
    
    def add_profile(self, profile: DoseResponseProfile):
        """
        Add a profile (should have S' calculated).
        
        Args:
            profile: DoseResponseProfile with S' calculated
            
        Raises:
            ValueError: If profile doesn't have S' or Hill params, or already exists
        """
        # Validation: profile should be processed
        if profile.s_prime is None:
            raise ValueError("Profile must have S' calculated before adding to ScreeningDataset")
        if profile.hill_params is None:
            raise ValueError("Profile must have Hill params before adding to ScreeningDataset")
        
        key = (profile.compound.drug_id, profile.cell_line.name)
        if key in self._profiles:
            raise ValueError(f"Profile for {key} already exists")
        self._profiles[key] = profile
    
    def get_profile(self, compound: Union[Compound, str], 
                    cell_line: Union[CellLine, str]) -> Optional[DoseResponseProfile]:
        """
        Retrieve a specific profile.
        
        Args:
            compound: Compound object or drug_id string
            cell_line: CellLine object or cell_line name string
            
        Returns:
            DoseResponseProfile or None if not found
        """
        compound_id = compound.drug_id if isinstance(compound, Compound) else compound
        cellline_name = cell_line.name if isinstance(cell_line, CellLine) else cell_line
        return self._profiles.get((compound_id, cellline_name))
    
    def calculate_delta_s_prime(
        self,
        reference_cell_lines: Union[str, List[str]],
        test_cell_lines: Union[str, List[str]],
        headings_one_to_one_in_ref_and_test: Optional[List[str]] = None,
        source_profile: Literal["ref", "test"] = "test",
    ) -> Dict[str, List[Dict]]:
        """
        Calculate delta S' = S'(ref) - S'(test) for each compound.
        
        Compound-level columns (1:1 per compound, not per cell line) auto-propagate:
        MOA and drug targets are reserved and always included; additional headings
        may be specified via headings_one_to_one_in_ref_and_test. Values are taken
        from the ref or test profile according to source_profile.
        
        Args:
            reference_cell_lines: Reference cell line name(s)
            test_cell_lines: Test cell line name(s)
            headings_one_to_one_in_ref_and_test: Optional list of metadata headings
                that exist 1:1 in ref and test; included in output, values from
                source_profile.
            source_profile: 'ref' or 'test'; which profile to use for compound-level
                values (MOA, drug targets, and optional headings).
            
        Returns:
            Dictionary with keys for each reference cell line, containing
            lists of dicts with delta S' and compound-level fields per combo.
        """
        ref_list = [reference_cell_lines] if isinstance(reference_cell_lines, str) else reference_cell_lines
        test_list = [test_cell_lines] if isinstance(test_cell_lines, str) else test_cell_lines
        extra_headings = headings_one_to_one_in_ref_and_test or []
        
        results = {}
        for ref_cellline in ref_list:
            rows = []
            compounds = {profile.compound.drug_id: profile.compound for profile in self._profiles.values()}
            
            for drug_id, compound in compounds.items():
                ref_profile = self.get_profile(compound, ref_cellline)
                if ref_profile is None or ref_profile.s_prime is None:
                    continue
                
                for test_cellline in test_list:
                    test_profile = self.get_profile(compound, test_cellline)
                    if test_profile is None or test_profile.s_prime is None:
                        continue
                    
                    delta = ref_profile.s_prime - test_profile.s_prime
                    source_meta = (ref_profile.metadata if source_profile == "ref" else test_profile.metadata) or {}
                    row = {
                        'compound_name': compound.name,
                        'drug_id': drug_id,
                        'reference_cell_line': ref_cellline,
                        'test_cell_line': test_cellline,
                        's_prime_ref': ref_profile.s_prime,
                        's_prime_test': test_profile.s_prime,
                        'delta_s_prime': delta,
                        'MOA': _resolve_moa(source_meta),
                        'drug targets': _resolve_drug_targets(source_meta),
                    }
                    for h in extra_headings:
                        row[h] = source_meta.get(h, "")
                    rows.append(row)
            
            results[ref_cellline] = rows
        
        return results
    
    def to_dict_list(self) -> List[Dict]:
        """
        Export to list of dictionaries with all S' values.
        
        Returns:
            List of dictionaries with profile data
        """
        rows = []
        for profile in self._profiles.values():
            row = {
                'compound_name': profile.compound.name,
                'drug_id': profile.compound.drug_id,
                'cell_line': profile.cell_line.name,
                's_prime': profile.s_prime,
                'ec50': profile.hill_params.ec50 if profile.hill_params else None,
                'upper': profile.hill_params.upper if profile.hill_params else None,
                'lower': profile.hill_params.lower if profile.hill_params else None,
                'rank': profile.rank,
            }
            if profile.hill_params:
                row['hill_coefficient'] = profile.hill_params.hill_coefficient
                row['r_squared'] = profile.hill_params.r_squared
            rows.append(row)
        return rows
    
    def export_to_csv(
        self,
        filepath: Union[str, Path],
        include_metadata: bool = True
    ) -> None:
        """
        Export all profiles to CSV file.
        
        Base columns (identifiers, Hill params, S', Rank) are always written.
        When include_metadata is True, all generic metadata keys (union across
        profiles) are included as pass-through columns.
        
        Args:
            filepath: Path to output CSV file
            include_metadata: When True, include all metadata columns. When False,
                only base columns.
        """
        filepath = Path(filepath)
        profiles = sorted(self.profiles, key=lambda p: (p.compound.name, p.cell_line.name))
        
        base_fieldnames = [
            'Compound Name', 'Compound_ID', 'pubchem_sid', 'SMILES',
            'Cell_Line', 'Cell_Line_Ref_ID',
            'EC50', 'Upper', 'Lower', 'Hill_Slope', 'r2',
            "S'", 'Rank'
        ]
        all_meta_keys = sorted(set(
            k for p in profiles if p.metadata for k in p.metadata
        ))
        fieldnames = list(base_fieldnames)
        if include_metadata and all_meta_keys:
            fieldnames.extend(all_meta_keys)
        
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for profile in profiles:
                meta = profile.metadata or {}
                row = {
                    'Compound Name': profile.compound.name,
                    'Compound_ID': profile.compound.drug_id,
                    'pubchem_sid': profile.compound.pubchem_sid or '',
                    'SMILES': profile.compound.smiles or '',
                    'Cell_Line': profile.cell_line.name,
                    'Cell_Line_Ref_ID': profile.cell_line.ref_id or '',
                    'EC50': f"{profile.hill_params.ec50:.6e}" if profile.hill_params else '',
                    'Upper': f"{profile.hill_params.upper:.2f}" if profile.hill_params else '',
                    'Lower': f"{profile.hill_params.lower:.2f}" if profile.hill_params else '',
                    'Hill_Slope': f"{profile.hill_params.hill_coefficient:.4f}" if profile.hill_params and profile.hill_params.hill_coefficient else '',
                    'r2': f"{profile.hill_params.r_squared:.4f}" if profile.hill_params and profile.hill_params.r_squared is not None else '',
                    "S'": f"{profile.s_prime:.4f}" if profile.s_prime else '',
                    'Rank': str(profile.rank) if profile.rank else '',
                }
                if include_metadata and all_meta_keys:
                    for k in all_meta_keys:
                        row[k] = meta.get(k, '')
                
                writer.writerow(row)
    
    @staticmethod
    def export_delta_s_prime_to_csv(
        delta_results: Dict[str, List[Dict]],
        filepath: Union[str, Path],
        headings_one_to_one_in_ref_and_test: Optional[List[str]] = None,
    ) -> None:
        """
        Export delta S' results to CSV file.
        
        Includes compound-level columns: MOA, drug targets (reserved), plus any
        headings specified in headings_one_to_one_in_ref_and_test. These must
        match the headings passed to calculate_delta_s_prime when producing
        delta_results.
        
        Args:
            delta_results: Dictionary from calculate_delta_s_prime()
            filepath: Path to output CSV file
            headings_one_to_one_in_ref_and_test: Optional list of metadata
                headings included in delta output (same as for calculate_delta_s_prime).
        """
        filepath = Path(filepath)
        extra_headings = headings_one_to_one_in_ref_and_test or []
        
        flat_results = []
        for ref_cellline, comparisons in delta_results.items():
            for comp in comparisons:
                row = {
                    'Compound Name': comp.get('compound_name', ''),
                    'Compound_ID': comp.get('drug_id', ''),
                    'Reference_Cell_Line': comp.get('reference_cell_line', ''),
                    'Test_Cell_Line': comp.get('test_cell_line', ''),
                    "S' (Reference)": f"{comp.get('s_prime_ref', 0.0):.4f}",
                    "S' (Test)": f"{comp.get('s_prime_test', 0.0):.4f}",
                    "Delta S'": f"{comp.get('delta_s_prime', 0.0):.4f}",
                    'MOA': comp.get('MOA', ''),
                    'drug targets': comp.get('drug targets', ''),
                }
                for h in extra_headings:
                    row[h] = comp.get(h, '')
                flat_results.append(row)
        
        flat_results.sort(key=lambda x: float(x["Delta S'"]))
        for rank, result in enumerate(flat_results, start=1):
            result['Rank'] = str(rank)
        
        fieldnames = [
            'Rank', 'Compound Name', 'Compound_ID',
            'Reference_Cell_Line', 'Test_Cell_Line',
            "S' (Reference)", "S' (Test)", "Delta S'",
            'MOA', 'drug targets',
        ]
        fieldnames.extend(extra_headings)
        
        with open(filepath, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(flat_results)
    
    @property
    def profiles(self):
        """Iterator over all profiles"""
        return self._profiles.values()
    
    def __len__(self):
        return len(self._profiles)


# ============================================================================
# Main API Class
# ============================================================================

class SPrime:
    """
    Main API for sprime package.
    
    Provides factory methods and convenience functions for loading data,
    processing dose-response profiles, and calculating S' values.
    """
    
    @staticmethod
    def load(
        filepath: Union[str, Path],
        assay_name: Optional[str] = None,
        values_as: str = "columns",
        **csv_kwargs
    ) -> Tuple[RawDataset, Optional['ProcessingReport']]:
        """
        Load raw data from CSV file or pandas DataFrame with quality reporting.
        
        Auto-detects pandas DataFrame input and converts it. For CSV files,
        uses global reporting configuration for console/log output.
        
        Args:
            filepath: Path to CSV file, or pandas DataFrame
            assay_name: Name for assay (defaults to filename stem or 'DataFrame')
            values_as: "columns" (DATA*/CONC*) or "list" (Responses, Concentrations)
            **csv_kwargs: Additional arguments for csv.DictReader (ignored for DataFrames)
            
        Returns:
            Tuple of (RawDataset, ProcessingReport)
            
        Raises:
            ValueError: If required columns are missing
            ImportError: If DataFrame provided but pandas not installed
        """
        # Auto-detect DataFrame
        try:
            import pandas as pd
            if isinstance(filepath, pd.DataFrame):
                return SPrime.load_from_dataframe(filepath, assay_name, values_as=values_as)
        except ImportError:
            pass  # pandas not available, treat as file path
        except (AttributeError, TypeError):
            pass  # Not a DataFrame, treat as file path
        
        # Treat as file path
        raw_dataset, report = RawDataset.load_from_file(
            filepath, assay_name, values_as=values_as, **csv_kwargs
        )
        
        # Auto-print and write log based on global config
        if report and ReportingConfig is not None:
            report.print_console_summary()
            report.write_log_file()
        
        return raw_dataset, report
    
    @staticmethod
    def load_from_dataframe(
        df,
        assay_name: Optional[str] = None,
        values_as: str = "columns",
    ) -> Tuple[RawDataset, Optional['ProcessingReport']]:
        """
        Load raw data from pandas DataFrame with quality reporting.
        
        Args:
            df: pandas DataFrame with columns matching CSV format
            assay_name: Name for assay (defaults to 'DataFrame')
            values_as: "columns" (DATA*/CONC*) or "list" (Responses, Concentrations)
            
        Returns:
            Tuple of (RawDataset, ProcessingReport)
            
        Raises:
            ValueError: If required columns are missing
            ImportError: If pandas is not installed
            TypeError: If input is not a pandas DataFrame
        """
        # Convert DataFrame to list of dicts
        list_of_rows = _convert_dataframe_to_dict_list(df)
        
        # Use get_s_prime_from_data logic but return RawDataset
        if not list_of_rows:
            # Create empty dataset
            if assay_name is None:
                assay_name = "DataFrame"
            assay = Assay(name=assay_name)
            raw_dataset = RawDataset(assay=assay)
            if ProcessingReport is not None:
                report = ProcessingReport()
            else:
                report = None
            return raw_dataset, report
        
        # Validate required columns exist
        first_row_keys = list(list_of_rows[0].keys()) if list_of_rows else []
        if first_row_keys:
            _validate_required_columns(
                first_row_keys, source_name="DataFrame", values_as=values_as
            )
        
        # Create report
        if ProcessingReport is not None:
            report = ProcessingReport()
        else:
            report = None
        
        # Create RawDataset
        if assay_name is None:
            assay_name = "DataFrame"
        assay = Assay(name=assay_name)
        raw_dataset = RawDataset(assay=assay)
        
        if report:
            report.total_rows = len(list_of_rows)
        
        # Process rows (similar to get_s_prime_from_data but add to RawDataset)
        # Rows are taken literally: empty values are null, no forward-filling.
        reserved = _reserved_column_names(values_as, first_row_keys)
        compounds_seen = set()
        
        for row_idx, row in enumerate(list_of_rows):
            if report:
                report.rows_processed += 1
            
            # Check if row is fully blank
            is_fully_blank = not any(
                v.strip() if isinstance(v, str) else str(v).strip()
                for v in row.values() if v
            )
            
            if is_fully_blank:
                continue
            
            # Check for empty cell line - RAISE EXCEPTION
            cell_line_name = row.get('Cell_Line', '').strip()
            if not cell_line_name:
                raise ValueError(
                    f"Row {row_idx + 1}: Missing required 'Cell_Line' value in DataFrame. "
                    f"All rows must have a cell line specified."
                )
            
            # Get compound info (Compound_ID required; NCGCID pass-through only)
            compound_name = row.get('Compound Name', '').strip() or 'Unknown'
            compound_id = row.get('Compound_ID', '').strip()
            
            if not compound_id:
                raise ValueError(
                    f"Row {row_idx + 1}: Missing required 'Compound_ID' value in DataFrame. "
                    f"Compound: {compound_name}, Cell_Line: {cell_line_name}. "
                    f"All rows must have a compound identifier."
                )
            
            # Track compound
            if report and compound_id not in compounds_seen:
                compounds_seen.add(compound_id)
                report.compounds_loaded += 1
            
            # Create compound and cell line objects
            compound = Compound(
                name=compound_name,
                drug_id=compound_id,
                pubchem_sid=row.get('pubchem_sid', '').strip() or None,
                smiles=row.get('SMILES', '').strip() or None
            )
            
            cell_line = CellLine(
                name=cell_line_name,
                ref_id=row.get('Cell_Line_Ref_ID', '').strip() or None
            )
            
            # Extract raw dose-response data (if present)
            concentrations = None
            responses = None
            
            if values_as == "list":
                resp_key = next((k for k in row if k.lower() == 'responses'), None)
                conc_key = next((k for k in row if k.lower() == 'concentrations'), None)
                if resp_key and conc_key:
                    resp_str = (row.get(resp_key) or '').strip()
                    conc_str = (row.get(conc_key) or '').strip()
                    if resp_str and conc_str:
                        responses = []
                        conc_parts = []
                        for part in resp_str.split(','):
                            t = part.strip()
                            if t:
                                try:
                                    v = float(t)
                                    if not (math.isnan(v) or math.isinf(v)):
                                        responses.append(v)
                                except (ValueError, TypeError):
                                    pass
                        for part in conc_str.split(','):
                            t = part.strip()
                            if t:
                                try:
                                    v = float(t)
                                    if not (math.isnan(v) or math.isinf(v)):
                                        conc_parts.append(v)
                                except (ValueError, TypeError):
                                    pass
                        if len(responses) == len(conc_parts) and len(responses) >= 4:
                            units = row.get('Concentration_Units', '').strip()
                            if not units:
                                raise ValueError(
                                    f"Row {row_idx + 1}: Missing required 'Concentration_Units' for raw data in DataFrame. "
                                    f"Compound: {compound_name}, Cell_Line: {cell_line_name}."
                                )
                            concentrations = convert_to_micromolar(conc_parts, units)
                        else:
                            if len(responses) != len(conc_parts):
                                raise ValueError(
                                    f"Row {row_idx + 1}: Responses and Concentrations length mismatch "
                                    f"({len(responses)} vs {len(conc_parts)}) in DataFrame."
                                )
                            responses = None
                            concentrations = None
            
            else:
                data_cols = [k for k in row.keys() if k.startswith('Data') or k.startswith('DATA')]
                conc_cols = [k for k in row.keys() 
                            if (k.startswith('Conc') or k.startswith('CONC')) 
                            and 'Units' not in k and 'units' not in k]
                if data_cols and conc_cols:
                    data_cols = sorted(data_cols, key=lambda x: int(''.join(filter(str.isdigit, x)) or '0'))
                    conc_cols = sorted(conc_cols, key=lambda x: int(''.join(filter(str.isdigit, x)) or '0'))
                    responses = []
                    concentrations = []
                    for data_col, conc_col in zip(data_cols, conc_cols):
                        try:
                            resp_val = row.get(data_col, '') or ''
                            conc_val = row.get(conc_col, '') or ''
                            resp_val = resp_val.strip() if isinstance(resp_val, str) else str(resp_val)
                            conc_val = conc_val.strip() if isinstance(conc_val, str) else str(conc_val)
                            if resp_val and conc_val:
                                responses.append(float(resp_val))
                                concentrations.append(float(conc_val))
                        except (ValueError, TypeError):
                            continue
                    if concentrations and responses:
                        units = row.get('Concentration_Units', '').strip()
                        if not units:
                            raise ValueError(
                                f"Row {row_idx + 1}: Missing required 'Concentration_Units' for raw data in DataFrame. "
                                f"Compound: {compound_name}, Cell_Line: {cell_line_name}. "
                                f"Raw dose-response data requires Concentration_Units."
                            )
                        concentrations = convert_to_micromolar(concentrations, units)
                    else:
                        concentrations = None
                        responses = None
            
            # Extract pre-calculated Hill params (if present)
            hill_params = None
            ac50 = row.get('AC50', '').strip() or row.get('ec50', '').strip()
            if ac50:
                try:
                    hill_params = HillCurveParams(
                        ec50=float(ac50),
                        upper=float(row.get('Upper', row.get('Infinity', '0')).strip() or '0'),
                        lower=float(row.get('Lower', row.get('Zero', '0')).strip() or '0'),
                        hill_coefficient=_try_float(row.get('Hill_Slope', row.get('Hill', row.get('slope', '')))),
                        r_squared=_try_float(row.get('r2', row.get('R²', '')))
                    )
                except (ValueError, TypeError):
                    hill_params = None
            
            # Validate that row has either raw data or pre-calculated params
            has_raw_data = concentrations is not None and responses is not None and len(concentrations) > 0
            has_precalc_params = hill_params is not None
            
            if not (has_raw_data or has_precalc_params):
                raise ValueError(
                    f"Row {row_idx + 1}: No dose-response data found for compound '{compound_name}' "
                    f"(Compound_ID: {compound_id}) in cell line '{cell_line_name}' in DataFrame. "
                    f"Row must have either: (1) raw data columns (DATA*/CONC*), or "
                    f"(2) pre-calculated parameters (AC50/Upper/Lower)."
                )
            
            # Pre-calculated S' (if present)
            s_prime = _try_float(row.get("S'", row.get('S Prime', '')))
            rank = _try_int(row.get('Rank', ''))
            
            # Extract metadata: generic pass-through (all non-reserved columns, exact header, value as-is)
            metadata = {}
            for col in row.keys():
                if col in reserved:
                    continue
                raw = row.get(col, "")
                raw = raw if isinstance(raw, str) else str(raw)
                metadata[col] = raw
            
            # Create profile
            profile = DoseResponseProfile(
                compound=compound,
                cell_line=cell_line,
                assay=assay,
                concentrations=concentrations,
                responses=responses,
                concentration_units="microM",
                hill_params=hill_params,
                s_prime=s_prime,
                rank=rank,
                metadata=metadata if metadata else None
            )
            
            raw_dataset.add_profile(profile)
            if report:
                report.profiles_created += 1
        
        # Auto-print and write log based on global config
        if report and ReportingConfig is not None:
            report.print_console_summary()
            report.write_log_file()
        
        return raw_dataset, report
    
    @staticmethod
    def process(
        raw_dataset: RawDataset,
        report: Optional['ProcessingReport'] = None,
        allow_overwrite_hill_coefficients: bool = False,
        **fit_params
    ) -> Tuple[ScreeningDataset, Optional['ProcessingReport']]:
        """
        Convert RawDataset to ScreeningDataset (fit curves, calculate S').
        
        Uses global reporting configuration for console/log output.
        
        Args:
            raw_dataset: RawDataset to process
            report: Optional ProcessingReport to accumulate warnings (reuses from load if None)
            allow_overwrite_hill_coefficients: If True, allow overwriting pre-calc
                Hill params (AC50, Upper, Lower, Hill_Slope, r2) with fitted values when
                both raw and pre-calc exist. Default False (raise). When True,
                overwrites are logged as warnings.
            **fit_params: Parameters for curve fitting (e.g. maxfev, bounds).
            
        Returns:
            Tuple of (ScreeningDataset, ProcessingReport)
        """
        screening_dataset, report = raw_dataset.to_screening_dataset(
            report=report,
            allow_overwrite_hill_coefficients=allow_overwrite_hill_coefficients,
            **fit_params
        )
        
        # Auto-print and write log based on global config
        if report and ReportingConfig is not None:
            report.print_console_summary()
            report.write_log_file()
        
        return screening_dataset, report


# ============================================================================
# Utility Functions (matching original pseudo-code)
# ============================================================================

def fit_hill_from_raw_data(
    response_array: List[float],
    concentration_array: List[float],
    concentration_units: str = "microM",
    **hill_calc_params
) -> HillCurveParams:
    """
    Fit Hill curve from raw response and concentration arrays.
    
    Matches original pseudo-code: sigleHillFromRawData()
    
    Args:
        response_array: List of response values
        concentration_array: List of concentration values
        concentration_units: Units of concentration (default: 'microM')
        **hill_calc_params: Additional parameters for curve fitting
        
    Returns:
        HillCurveParams: Fitted curve parameters
    """
    if hill_fitting is None:
        raise ImportError(
            "Hill curve fitting requires scipy. "
            "Install with: pip install scipy"
        )
    
    # Convert units if needed
    concentrations = convert_to_micromolar(concentration_array, concentration_units)
    
    # Fit using hill_fitting module
    return hill_fitting.fit_hill_curve(
        concentrations,
        response_array,
        **hill_calc_params
    )


def calculate_s_prime_from_params(ac50: float, upper: float, lower: float) -> float:
    """
    Calculate S' from Hill curve parameters.
    
    Matches original pseudo-code: singleSPrime()
    
    Args:
        ac50: AC50 or EC50 value
        upper: Upper asymptote
        lower: Lower asymptote
        
    Returns:
        float: S' value
    """
    # S' = asinh((Upper-Lower)/EC50)
    ratio = (upper - lower) / ac50
    return math.asinh(ratio)


def get_s_primes_from_file(
    filepath: Union[str, Path],
    allow_overwrite_hill_coefficients: bool = False,
    values_as: str = "columns",
    **fit_params
) -> List[Dict]:
    """
    Load CSV and calculate S' values for all compounds.
    
    Matches original pseudo-code: getSPrimesFromFile()
    Uses global reporting configuration for console/log output.
    
    Args:
        filepath: Path to CSV file
        allow_overwrite_hill_coefficients: If True, allow overwriting pre-calc
            Hill params (AC50, Upper, Lower, Hill, r2) with fitted values when
            both raw and pre-calc exist. Default False (raise). When True,
            overwrites are logged as warnings.
        values_as: "columns" (DATA*/CONC*) or "list" (Responses, Concentrations)
        **fit_params: Parameters for curve fitting (e.g. maxfev, bounds).
        
    Returns:
        List of dictionaries with S' values and ranking
    """
    # Create single report that accumulates both load and process warnings
    if ProcessingReport is not None:
        report = ProcessingReport()
        report.input_filepath = Path(filepath)
    else:
        report = None
    
    # Load
    raw_dataset, load_report = RawDataset.load_from_file(
        filepath, report=report, values_as=values_as
    )
    
    # Process (reuses same report)
    screening_dataset, process_report = raw_dataset.to_screening_dataset(
        report=report,
        allow_overwrite_hill_coefficients=allow_overwrite_hill_coefficients,
        **fit_params
    )
    
    # Print and write log based on global config
    if report and ReportingConfig is not None:
        report.print_console_summary()
        report.write_log_file()
    
    # Add ranking
    profiles = list(screening_dataset.profiles)
    profiles.sort(key=lambda p: p.s_prime if p.s_prime is not None else float('-inf'), reverse=True)
    
    for rank, profile in enumerate(profiles, start=1):
        profile.rank = rank
    
    return screening_dataset.to_dict_list()


def get_s_prime_from_data(
    list_of_rows: List[Dict],
    allow_overwrite_hill_coefficients: bool = False,
    values_as: str = "columns",
    **fit_params
) -> List[Dict]:
    """
    Calculate S' values from in-memory data structure.
    
    Matches original pseudo-code: getSPrimeFromData()
    Uses global reporting configuration for console/log output.
    
    Args:
        list_of_rows: List of dictionaries matching CSV row format
        allow_overwrite_hill_coefficients: If True, allow overwriting pre-calc
            Hill params (AC50, Upper, Lower, Hill, r2) with fitted values when
            both raw and pre-calc exist. Default False (raise). When True,
            overwrites are logged as warnings.
        values_as: "columns" (DATA*/CONC*) or "list" (Responses, Concentrations)
        **fit_params: Parameters for curve fitting (e.g. maxfev, bounds).
        
    Returns:
        List of dictionaries with S' values and ranking
        
    Raises:
        ValueError: If required columns are missing or data is invalid
    """
    if not list_of_rows:
        return []
    
    # Validate required columns exist (check first row keys)
    first_row_keys = list(list_of_rows[0].keys()) if list_of_rows else []
    if first_row_keys:
        _validate_required_columns(
            first_row_keys, source_name="in-memory data", values_as=values_as
        )
    
    # Create report (row numbers will be 0 for in-memory data)
    if ProcessingReport is not None:
        report = ProcessingReport()
    else:
        report = None
    
    # Create a temporary RawDataset by manually constructing profiles
    # We'll use the first row to infer assay name, or use a default
    assay_name = list_of_rows[0].get('Screen ID', list_of_rows[0].get('Assay', 'Unknown'))
    assay = Assay(name=assay_name)
    raw_dataset = RawDataset(assay=assay)
    
    if report:
        report.total_rows = len(list_of_rows)
    
    reserved = _reserved_column_names(values_as, first_row_keys)
    compounds_seen = set()
    
    for row_idx, row in enumerate(list_of_rows):
        if report:
            report.rows_processed += 1
        
        # Check if row is fully blank
        is_fully_blank = not any(
            v.strip() if isinstance(v, str) else str(v).strip()
            for v in row.values() if v
        )
        
        if is_fully_blank:
            continue
        
        # Get cell line name early for use in warnings
        cell_line_name = row.get('Cell_Line', '').strip()
        
        # Check for empty cell line - RAISE EXCEPTION
        if not cell_line_name:
            raise ValueError(
                f"Row {row_idx + 1}: Missing required 'Cell_Line' value in in-memory data. "
                f"All rows must have a cell line specified."
            )
        
        # Get compound info (Compound_ID required; NCGCID pass-through only)
        compound_name = row.get('Compound Name', '').strip() or 'Unknown'
        compound_id = row.get('Compound_ID', '').strip()
        
        if not compound_id:
            raise ValueError(
                f"Row {row_idx + 1}: Missing required 'Compound_ID' value in in-memory data. "
                f"Compound: {compound_name}, Cell_Line: {cell_line_name}. "
                f"All rows must have a compound identifier."
            )
        
        # Track compound
        if report and compound_id not in compounds_seen:
            compounds_seen.add(compound_id)
            report.compounds_loaded += 1
        
        # Check for missing compound name
        if report and compound_name == 'Unknown':
            report.add_warning(
                row_number=0,
                category="DATA_QUALITY",
                message="Compound Name missing, using 'Unknown'",
                drug_id=compound_id,
                cell_line=cell_line_name,
                field_name="Compound Name"
            )
            report.missing_compound_names += 1
        
        compound = Compound(
            name=compound_name,
            drug_id=compound_id,
            pubchem_sid=row.get('pubchem_sid', '').strip() or None,
            smiles=row.get('SMILES', '').strip() or None
        )
        
        # Create CellLine (cell_line_name already defined above)
        cell_line = CellLine(
            name=cell_line_name,
            ref_id=row.get('Cell_Line_Ref_ID', '').strip() or None
        )
        
        # Extract raw dose-response data (if present)
        concentrations = None
        responses = None
        
        if values_as == "list":
            resp_key = next((k for k in row if k.lower() == 'responses'), None)
            conc_key = next((k for k in row if k.lower() == 'concentrations'), None)
            if resp_key and conc_key:
                resp_str = (row.get(resp_key) or '').strip()
                conc_str = (row.get(conc_key) or '').strip()
                if resp_str and conc_str:
                    responses = []
                    conc_parts = []
                    for part in resp_str.split(','):
                        t = part.strip()
                        if t:
                            try:
                                v = float(t)
                                if not (math.isnan(v) or math.isinf(v)):
                                    responses.append(v)
                            except (ValueError, TypeError):
                                pass
                    for part in conc_str.split(','):
                        t = part.strip()
                        if t:
                            try:
                                v = float(t)
                                if not (math.isnan(v) or math.isinf(v)):
                                    conc_parts.append(v)
                            except (ValueError, TypeError):
                                pass
                    if len(responses) == len(conc_parts) and len(responses) >= 4:
                        units = row.get('Concentration_Units', '').strip()
                        if not units:
                            raise ValueError(
                                f"Row {row_idx + 1}: Missing required 'Concentration_Units' for raw data in in-memory data. "
                                f"Compound: {compound_name}, Cell_Line: {cell_line_name}."
                            )
                        concentrations = convert_to_micromolar(conc_parts, units)
                    else:
                        if len(responses) != len(conc_parts):
                            raise ValueError(
                                f"Row {row_idx + 1}: Responses and Concentrations length mismatch "
                                f"({len(responses)} vs {len(conc_parts)}) in in-memory data."
                            )
                        responses = None
                        concentrations = None
        
        else:
            data_cols = [k for k in row.keys() if k.startswith('Data') or k.startswith('DATA')]
            conc_cols = [k for k in row.keys() 
                        if (k.startswith('Conc') or k.startswith('CONC')) 
                        and 'Units' not in k and 'units' not in k]
            if data_cols and conc_cols:
                data_cols = sorted(data_cols, key=lambda x: int(''.join(filter(str.isdigit, x)) or '0'))
                conc_cols = sorted(conc_cols, key=lambda x: int(''.join(filter(str.isdigit, x)) or '0'))
                responses = []
                concentrations = []
                for data_col, conc_col in zip(data_cols, conc_cols):
                    try:
                        resp_val = row.get(data_col, '') or ''
                        conc_val = row.get(conc_col, '') or ''
                        resp_val = resp_val.strip() if isinstance(resp_val, str) else ''
                        conc_val = conc_val.strip() if isinstance(conc_val, str) else ''
                        if resp_val and conc_val:
                            responses.append(float(resp_val))
                            concentrations.append(float(conc_val))
                    except (ValueError, TypeError):
                        continue
                if concentrations and responses:
                    units = row.get('Concentration_Units', '').strip()
                    if not units:
                        raise ValueError(
                            f"Row {row_idx + 1}: Missing required 'Concentration_Units' for raw data in in-memory data. "
                            f"Compound: {compound_name}, Cell_Line: {cell_line_name}."
                        )
                    concentrations = convert_to_micromolar(concentrations, units)
                else:
                    concentrations = None
                    responses = None
        
        # Extract pre-calculated Hill params (if present)
        hill_params = None
        ac50 = row.get('AC50', '').strip() or row.get('ec50', '').strip()
        if ac50:
            try:
                hill_params = HillCurveParams(
                    ec50=float(ac50),
                    upper=float(row.get('Upper', row.get('Infinity', '0')).strip() or '0'),
                    lower=float(row.get('Lower', row.get('Zero', '0')).strip() or '0'),
                    hill_coefficient=_try_float(row.get('Hill_Slope', row.get('Hill', row.get('slope', '')))),
                    r_squared=_try_float(row.get('r2', row.get('R²', '')))
                )
            except (ValueError, TypeError):
                hill_params = None
        
        # Validate that row has either raw data or pre-calculated params
        has_raw_data = concentrations is not None and responses is not None and len(concentrations) > 0
        has_precalc_params = hill_params is not None
        
        if not (has_raw_data or has_precalc_params):
            raise ValueError(
                f"Row {row_idx + 1}: No dose-response data found for compound '{compound_name}' "
                f"(Compound_ID: {compound_id}) in cell line '{cell_line_name}' in in-memory data. "
                f"Row must have either: (1) raw data columns (DATA*/CONC*), or "
                f"(2) pre-calculated parameters (AC50/Upper/Lower)."
            )
        
        # Pre-calculated S' (if present)
        s_prime = _try_float(row.get("S'", row.get('S Prime', '')))
        rank = _try_int(row.get('Rank', ''))
        
        # Extract metadata: generic pass-through (all non-reserved columns, exact header, value as-is)
        metadata = {}
        for col in row.keys():
            if col in reserved:
                continue
            raw = row.get(col, "")
            raw = raw if isinstance(raw, str) else str(raw)
            metadata[col] = raw
        
        # Create profile
        profile = DoseResponseProfile(
            compound=compound,
            cell_line=cell_line,
            assay=assay,
            concentrations=concentrations,
            responses=responses,
            concentration_units="microM",
            hill_params=hill_params,
            s_prime=s_prime,
            rank=rank,
            metadata=metadata if metadata else None
        )
        
        raw_dataset.add_profile(profile)
        if report:
            report.profiles_created += 1
    
    # Process to ScreeningDataset (reuses same report)
    screening_dataset, process_report = raw_dataset.to_screening_dataset(
        report=report,
        allow_overwrite_hill_coefficients=allow_overwrite_hill_coefficients,
        **fit_params
    )
    
    # Print and write log based on global config
    if report and ReportingConfig is not None:
        report.print_console_summary()
        report.write_log_file()
    
    # Add ranking
    profiles = list(screening_dataset.profiles)
    profiles.sort(key=lambda p: p.s_prime if p.s_prime is not None else float('-inf'), reverse=True)
    
    for rank, profile in enumerate(profiles, start=1):
        profile.rank = rank
    
    return screening_dataset.to_dict_list()


def calculate_delta_s_prime(
    s_prime_data: Union[ScreeningDataset, List[Dict]],
    reference_cell_line_names: Union[str, List[str]],
    test_cell_line_names: Union[str, List[str]],
    headings_one_to_one_in_ref_and_test: Optional[List[str]] = None,
    source_profile: Literal["ref", "test"] = "test",
) -> Dict[str, List[Dict]]:
    """
    Calculate delta S' between reference and test cell lines.
    
    Matches original pseudo-code: delta_s_prime()
    
    Compound-level columns (MOA, drug targets, optional headings) auto-propagate;
    see ScreeningDataset.calculate_delta_s_prime for details.
    
    Args:
        s_prime_data: ScreeningDataset or list of dicts with S' values
        reference_cell_line_names: Reference cell line name(s)
        test_cell_line_names: Test cell line name(s)
        headings_one_to_one_in_ref_and_test: Optional list of metadata headings
            that exist 1:1 in ref and test; included in output.
        source_profile: 'ref' or 'test'; which profile to use for compound-level values.
        
    Returns:
        Dictionary with keys for each reference cell line, containing
        lists of dicts with delta S' and compound-level fields per combo.
    """
    if isinstance(s_prime_data, ScreeningDataset):
        return s_prime_data.calculate_delta_s_prime(
            reference_cell_line_names,
            test_cell_line_names,
            headings_one_to_one_in_ref_and_test=headings_one_to_one_in_ref_and_test,
            source_profile=source_profile,
        )
    else:
        if not s_prime_data:
            return {}
        
        assay_name = s_prime_data[0].get('assay', 'Unknown')
        assay = Assay(name=assay_name)
        screening_dataset = ScreeningDataset(assay=assay)
        
        for row in s_prime_data:
            compound = Compound(
                name=row.get('compound_name', 'Unknown'),
                drug_id=row.get('drug_id', ''),
                pubchem_sid=None,
                smiles=None
            )
            cell_line = CellLine(
                name=row.get('cell_line', ''),
                ref_id=None
            )
            
            hill_params = None
            if row.get('ec50') is not None:
                hill_params = HillCurveParams(
                    ec50=row.get('ec50'),
                    upper=row.get('upper'),
                    lower=row.get('lower'),
                    hill_coefficient=row.get('hill_coefficient'),
                    r_squared=row.get('r_squared')
                )
            
            profile = DoseResponseProfile(
                compound=compound,
                cell_line=cell_line,
                assay=assay,
                concentrations=None,
                responses=None,
                concentration_units="microM",
                hill_params=hill_params,
                s_prime=row.get('s_prime'),
                rank=row.get('rank'),
                metadata=None
            )
            
            if profile.s_prime is not None and profile.hill_params is not None:
                try:
                    screening_dataset.add_profile(profile)
                except ValueError:
                    continue
        
        return screening_dataset.calculate_delta_s_prime(
            reference_cell_line_names,
            test_cell_line_names,
            headings_one_to_one_in_ref_and_test=headings_one_to_one_in_ref_and_test,
            source_profile=source_profile,
        )


# ============================================================================
# Helper Functions
# ============================================================================

def convert_to_micromolar(concentrations: List[float], units: str) -> List[float]:
    """
    Convert concentration values to microMolar.
    
    Args:
        concentrations: List of concentration values
        units: Current units. Supported (case-insensitive), smallest to largest:
            fM (fm, femtom); pM (pm, picom); nM (nm, nanom);
            microM (µm, um, microm, micro); mM (mm, millim); M (m, mol).
        
    Returns:
        List of concentrations in microMolar
    """
    units_lower = units.lower().strip()
    
    conversion_factors = {
        'fm': 1e-9,
        'femtom': 1e-9,
        'pm': 1e-6,
        'picom': 1e-6,
        'nm': 0.001,
        'nanom': 0.001,
        'microm': 1.0,
        'micro': 1.0,
        'μm': 1.0,
        'um': 1.0,
        'mm': 1000.0,
        'millim': 1000.0,
        'm': 1000000.0,
        'mol': 1000000.0,
    }
    
    factor = conversion_factors.get(units_lower, 1.0)
    return [c * factor for c in concentrations]


def _try_float(value: str) -> Optional[float]:
    """Try to convert string to float, return None if fails"""
    if not value or not value.strip():
        return None
    try:
        return float(value.strip())
    except (ValueError, TypeError):
        return None


def _try_int(value: str) -> Optional[int]:
    """Try to convert string to int, return None if fails"""
    if not value or not value.strip():
        return None
    try:
        return int(float(value.strip()))  # Convert via float to handle "3.0" -> 3
    except (ValueError, TypeError):
        return None
