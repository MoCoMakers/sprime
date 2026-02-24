"""
sprime - A biomedical library for screening high-throughput screening data.

sprime provides tools for analyzing and processing high-throughput screening
data in preclinical drug studies.
"""

try:
    from ._version import __version__
except ImportError:
    # Fallback for when version file doesn't exist (development without build)
    try:
        from importlib.metadata import version
        __version__ = version("sprime")
    except Exception:
        __version__ = "0.0.0.dev0+unknown"

__author__ = "MoCo Makers"

from .sprime import (
    # Core classes
    SPrime,
    RawDataset,
    ScreeningDataset,
    DoseResponseProfile,
    # Value objects
    Compound,
    CellLine,
    Assay,
    HillCurveParams,
    # Utility functions
    fit_hill_from_raw_data,
    calculate_s_prime_from_params,
    get_s_primes_from_file,
    get_s_prime_from_data,
    calculate_delta_s_prime,
    convert_to_micromolar,
    normalize_responses,
)
from . import hill_fitting

# Alias for alternative capitalization
Sprime = SPrime

# Import reporting classes
try:
    from .reporting import ReportingConfig, ConsoleOutput, ProcessingReport
except ImportError:
    ReportingConfig = None
    ConsoleOutput = None
    ProcessingReport = None

__all__ = [
    "SPrime", "Sprime",  # Both capitalizations work
    "RawDataset", "ScreeningDataset", "DoseResponseProfile",
    "Compound", "CellLine", "Assay", "HillCurveParams",
    "fit_hill_from_raw_data", "calculate_s_prime_from_params",
    "get_s_primes_from_file", "get_s_prime_from_data", "calculate_delta_s_prime",
    "convert_to_micromolar", "normalize_responses",
    "hill_fitting",
]

# Add reporting to exports if available
if ReportingConfig is not None:
    __all__.extend(["ReportingConfig", "ConsoleOutput", "ProcessingReport"])

