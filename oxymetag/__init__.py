"""
OxyMetaG: Oxygen metabolism profiling from metagenomic data
"""

__version__ = "1.1.2"
__author__ = "Clifton P. Bueno de Mesquita"
__email__ = "cliff.buenodemesquita@colorado.edu"

from .core import extract_reads, profile_samples, predict_aerobes
from .utils import check_dependencies, run_kraken2_setup

__all__ = [
    "extract_reads",
    "profile_samples", 
    "predict_aerobes",
    "check_dependencies",
    "run_kraken2_setup"
]
