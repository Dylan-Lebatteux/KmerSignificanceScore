"""
K-mer Significance Score (KSS) Package

A comprehensive toolkit for computing k-mer significance scores in genomic sequences,
combining mutational impact, discriminative power, and protein functional importance.

Main modules:
    - kss: Core KSS computation and compilation
    - mutation_score: Amino acid substitution scoring
    - discriminative_score: Class discrimination metrics
    - protein_score: Protein functional importance from UniProt
    - kanalyzer: Sequence alignment and mutation identification
    - utils: Data I/O utilities
"""

__version__ = "2.0.0"
__author__ = "KSS Development Team"

# Import main functions for convenient top-level access
from .kss import compute_kss_scores, compile_results, get_taxon_id
from .kanalyzer import analyze_records
from .utils import load_data_from_json, save_data_as_json

__all__ = [
    "compute_kss_scores",
    "compile_results",
    "get_taxon_id",
    "analyze_records",
    "load_data_from_json",
    "save_data_as_json",
]
