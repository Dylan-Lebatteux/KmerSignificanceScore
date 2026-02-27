#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Amino Acid Substitution Scoring Module

This module provides functions for computing mutational impact scores based on
amino acid substitution matrices. It evaluates the functional impact of amino
acid changes by leveraging biochemical distance metrics.

The module supports multiple substitution matrices including:
- MIYATA: Evolutionary distance based on physicochemical properties
- BLOSUM: Block substitution matrices from protein alignments
- GRANTHAM: Chemical distance between amino acids
- PAM: Point accepted mutation matrices

Features:
- Automatic format detection (JSON/Pickle)
- Matrix normalization and scaling to [0, 1] range
- Support for custom substitution matrices
- Handling of deletions and insertions
- High-precision value preservation

Main functions:
    get_mutational_scores: Calculate impact scores for amino acid changes
    load_substitution_matrix: Load matrix from file with auto-detection
    save_substitution_matrix: Save matrix in multiple formats
    adjust_and_scale_substitution_matrix: Normalize matrix values
"""

import os
import json
import pickle
import numpy as np
from typing import Dict, Any, List, Tuple, Optional, Union


# Cache for substitution matrices to avoid repeated file I/O
_MATRIX_CACHE = {}


def load_substitution_matrix(filename: str) -> Dict[Tuple[str, str], float]:
    """
    Load a substitution matrix from JSON or Pickle file with automatic format detection.

    The function first attempts to determine the file format from the extension,
    then falls back to trying both formats if extension is ambiguous.

    Args:
        filename: Path to the matrix file (.json or .pkl extension preferred)

    Returns:
        Substitution matrix as a dictionary with amino acid pair tuples as keys
        and substitution scores as values. Example: {('A', 'G'): 0.45, ...}

    Raises:
        FileNotFoundError: If the specified file doesn't exist
        RuntimeError: If the file format is not supported or corrupted

    Notes:
        - Pickle format is preferred for higher numerical precision
        - JSON keys are automatically converted from strings to tuples
        - Supports both ('A', 'G') and 'AG' key formats in JSON
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    
    # Automatic format detection by extension
    _, ext = os.path.splitext(filename.lower())

    if ext == '.pkl':
        return _load_matrix_pickle(filename)
    elif ext == '.json':
        return _load_matrix_json(filename)
    else:
        # Try automatic detection if extension is unknown
        try:
            return _load_matrix_pickle(filename)
        except (pickle.UnpicklingError, EOFError, AttributeError, ImportError):
            # Pickle format failed, try JSON
            return _load_matrix_json(filename)


def _load_matrix_pickle(filename: str) -> Dict[Tuple[str, str], float]:
    """
    Load substitution matrix from pickle file.

    Args:
        filename: Path to pickle (.pkl) file

    Returns:
        Substitution matrix dictionary

    Raises:
        RuntimeError: If pickle loading fails
    """
    try:
        with open(filename, 'rb') as f:
            matrix = pickle.load(f)
        return matrix
    except Exception as e:
        raise RuntimeError(f"Error loading pickle file {filename}: {e}")


def _load_matrix_json(filename: str) -> Dict[Tuple[str, str], float]:
    """
    Load substitution matrix from JSON file.

    Args:
        filename: Path to JSON (.json) file

    Returns:
        Substitution matrix dictionary with tuple keys

    Raises:
        RuntimeError: If JSON loading or parsing fails

    Notes:
        - Converts string keys to tuple format
        - Handles both '(A, G)' and 'AG' key formats
    """
    try:
        with open(filename, 'r') as f:
            matrix_data = json.load(f)
        
        matrix = {}
        for key, value in matrix_data.items():
            if ',' in key:
                aa_pair = tuple(aa.strip("'\" ()") for aa in key.strip().split(','))
            else:
                aa_pair = tuple(key.strip())
            
            matrix[aa_pair] = float(value)
        
        return matrix
        
    except (FileNotFoundError, json.JSONDecodeError) as e:
        raise RuntimeError(f"Error loading JSON file {filename}: {e}")


def save_substitution_matrix(matrix: Dict[Tuple[str, str], float],
                            base_filename: str,
                            formats: List[str] = ['pickle']) -> Dict[str, str]:
    """
    Save a substitution matrix in one or multiple formats.

    Args:
        matrix: Substitution matrix dictionary with tuple keys
        base_filename: Base filename without extension (e.g., 'matrices/BLOSUM62')
        formats: List of output formats. Options: 'pickle', 'json', or 'both'
                Default: ['pickle']

    Returns:
        Dictionary mapping format names to created file paths.
        Example: {'pickle': 'matrices/BLOSUM62.pkl', 'json': 'matrices/BLOSUM62.json'}

    Notes:
        - Creates parent directories if they don't exist
        - 'both' is converted to ['pickle', 'json']
        - JSON format uses 15 decimal places for precision
        - Pickle format preserves full Python float precision
    """
    if 'both' in formats:
        formats = ['pickle', 'json']
    
    created_files = {}
    
    # Create directory if necessary
    os.makedirs(os.path.dirname(base_filename), exist_ok=True)
    
    if 'pickle' in formats:
        pickle_path = f"{base_filename}.pkl"
        _save_matrix_pickle(matrix, pickle_path)
        created_files['pickle'] = pickle_path
    
    if 'json' in formats:
        json_path = f"{base_filename}.json"
        _save_matrix_json(matrix, json_path)
        created_files['json'] = json_path
    
    return created_files


def _save_matrix_pickle(matrix: Dict[Tuple[str, str], float], filepath: str) -> None:
    """
    Save matrix in pickle format.

    Args:
        matrix: Substitution matrix dictionary
        filepath: Complete path including .pkl extension
    """
    with open(filepath, 'wb') as f:
        pickle.dump(matrix, f)


def _save_matrix_json(matrix: Dict[Tuple[str, str], float], filepath: str) -> None:
    """
    Save matrix in JSON format with string keys.

    Args:
        matrix: Substitution matrix dictionary
        filepath: Complete path including .json extension

    Notes:
        - Converts tuple keys to string format '(A, G)'
        - Rounds values to 15 decimal places
    """
    string_matrix = {}
    for (aa1, aa2), value in matrix.items():
        key = f"({aa1}, {aa2})"
        string_matrix[key] = round(float(value), 15)
    
    with open(filepath, "w") as outfile:
        json.dump(string_matrix, outfile, indent=4, sort_keys=True)


def adjust_and_scale_substitution_matrix(matrix: Dict[Tuple[str, str], float],
                                        matrix_type: str) -> Dict[Tuple[str, str], float]:
    """
    Adjust and scale a substitution matrix to [0, 1] range.

    Distance-based matrices (MIYATA, GRANTHAM, etc.) are inverted so that
    similar amino acids receive higher scores. All matrices are then scaled
    using min-max normalization to ensure values fall in [0, 1] range.

    Args:
        matrix: Original substitution matrix with tuple keys
        matrix_type: Name/type of the matrix (e.g., "MIYATA", "BLOSUM62")
                    Used to determine if inversion is needed

    Returns:
        Scaled substitution matrix with values in [0, 1] range.
        Higher values indicate more conservative substitutions (or smaller
        distances for inverted matrices).

    Notes:
        - Distance matrices are inverted: new_value = max_value - old_value
        - Ensures symmetry: if (A,G) exists, (G,A) will have same value
        - Uses global min-max normalization to preserve symmetry
        - Distance matrix types: MIYATA, GRANTHAM, SNEATH, and variants
    """
    # Handle distance-based matrices by inverting values
    distance_matrices = {
        "MIYATA", "GRANTHAM", "SNEATH", "MIYATA_EVO", 
        "MIYATA_BEST_GLOBAL", "MIYATA_OPTIMIZED"
    }
    
    if matrix_type.upper() in distance_matrices:
        max_value = max(matrix.values())
        matrix = {key: max_value - value for key, value in matrix.items()}
    
    # Extract all amino acids from the matrix
    amino_acids = sorted(set(aa for pair in matrix for aa in pair))
    
    # Build a 2D numpy array for scaling
    matrix_array = np.zeros((len(amino_acids), len(amino_acids)))
    aa_to_index = {aa: i for i, aa in enumerate(amino_acids)}
    
    # Fill the array with substitution values
    for (aa1, aa2), value in matrix.items():
        i, j = aa_to_index[aa1], aa_to_index[aa2]
        matrix_array[i, j] = value
        # Ensure symmetry if not already present
        matrix_array[j, i] = value
    
    # Scale values to [0, 1] range using global min-max normalization
    # Global scaling preserves symmetry: score(A→G) == score(G→A)
    min_val = matrix_array.min()
    max_val = matrix_array.max()
    scaled_array = (matrix_array - min_val) / (max_val - min_val)
    
    # Convert back to dictionary
    scaled_matrix = {}
    for aa1 in amino_acids:
        for aa2 in amino_acids:
            i, j = aa_to_index[aa1], aa_to_index[aa2]
            scaled_matrix[(aa1, aa2)] = scaled_array[i, j]
    
    return scaled_matrix


def get_substitution_matrix(matrix_type: str) -> Dict[Tuple[str, str], float]:
    """
    Get a processed substitution matrix with automatic format selection.

    Attempts to load the matrix from the substitution_matrices directory,
    preferring pickle format over JSON to preserve numerical precision.

    Args:
        matrix_type: Name of the matrix to load (e.g., 'MIYATA_EVO', 'BLOSUM62')

    Returns:
        Processed substitution matrix with values scaled to [0, 1] range

    Raises:
        FileNotFoundError: If no matrix file (.pkl or .json) is found for the
                          specified matrix_type

    Notes:
        - Priority order: .pkl > .json (to avoid precision loss)
        - Automatically applies adjustment and scaling
        - Looks for files in the 'substitution_matrices' subdirectory
          relative to this module's location
        - Results are cached to avoid repeated file I/O
    """
    # OPTIMIZATION: Check cache first to avoid repeated file loading
    if matrix_type in _MATRIX_CACHE:
        return _MATRIX_CACHE[matrix_type]

    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Priority 1: Try pickle first (higher precision)
    pickle_path = os.path.join(current_dir, f'substitution_matrices/{matrix_type}.pkl')
    if os.path.exists(pickle_path):
        raw_matrix = load_substitution_matrix(pickle_path)
    else:
        # Priority 2: Fallback to JSON
        json_path = os.path.join(current_dir, f'substitution_matrices/{matrix_type}.json')
        if os.path.exists(json_path):
            raw_matrix = load_substitution_matrix(json_path)
        else:
            raise FileNotFoundError(f"No file found for {matrix_type} (.pkl or .json)")

    processed_matrix = adjust_and_scale_substitution_matrix(raw_matrix, matrix_type)

    # Cache the processed matrix for future use
    _MATRIX_CACHE[matrix_type] = processed_matrix

    return processed_matrix


def get_mutational_scores(changes: Union[List[str], Dict[str, Any]],
                         substitution_matrix_type: str = "MIYATA_EVO",
                         custom_matrix: Optional[Dict[Tuple[str, str], float]] = None,
                         min_score: float = 0.1,
                         indel_score: float = 1.0) -> Dict[str, float]:
    """
    Compute mutational impact scores for amino acid changes.

    This function evaluates the functional impact of amino acid mutations using
    substitution matrices. The scores reflect the biochemical dissimilarity
    between original and mutated amino acids.

    Args:
        changes: Amino acid changes to score. Can be:
                - List of mutation strings (e.g., ["A123G", "R45K", "D100-"])
                - Dictionary with mutation strings as keys
                Format: {original_aa}{position}{mutated_aa}
                Deletions use '-' as the mutated amino acid
        substitution_matrix_type: Name of substitution matrix to use.
                                 Default: "MIYATA_EVO"
                                 Ignored if custom_matrix is provided
                                 Options: MIYATA_EVO, BLOSUM62, GRANTHAM, etc.
        custom_matrix: Optional custom substitution matrix to use instead of
                      loading from file. Should have amino acid pair tuples as keys
        min_score: Minimum score for any substitution. Ensures even conservative
                  mutations have measurable impact. Default: 0.1
        indel_score: Score assigned to insertions and deletions (indels).
                    Default: 1.0 (maximum impact). Can be adjusted to
                    modulate the penalty for indel events.

    Returns:
        Dictionary mapping mutation strings to impact scores.
        Score range: [min_score, indel_score]
        - min_score: Most conservative/similar substitution
        - indel_score: Score for deletion/insertion events

    Notes:
        - Deletions ('-' in mutation string) receive the indel_score value
        - Mutation format: first char = original AA, last char = mutated AA
        - Uses inverse of similarity: score = 1.0 - similarity
        - Symmetric lookup: (A,G) and (G,A) have same score
        - Scores are rounded to 3 decimal places
    """
    # Handle both list and dict inputs
    if isinstance(changes, dict):
        mutation_list = list(changes.keys())
    else:
        mutation_list = changes
    
    # Get processed substitution matrix
    if custom_matrix is not None:
        # Use the provided custom matrix and process it
        matrix = adjust_and_scale_substitution_matrix(custom_matrix, substitution_matrix_type)
    else:
        # Load and process matrix from file (priority: pickle > json)
        matrix = get_substitution_matrix(substitution_matrix_type)

    # Calculate impact scores for each mutation
    scores = {}
    for mutation in mutation_list:
        # Assign indel_score for deletions/insertions
        if "-" in mutation:
            score = indel_score
        else:
            # Extract original and mutated amino acids
            original_aa = mutation[0]
            mutated_aa = mutation[-1]
            
            # Get substitution score (symmetric lookup)
            pair = (original_aa, mutated_aa)
            reverse_pair = (mutated_aa, original_aa)
            
            substitution_score = matrix.get(pair, matrix.get(reverse_pair, 0))
            
            # Calculate mutation score (higher = more significant change)
            raw_score = 1.0 - substitution_score
            
            # Apply minimum score to avoid zero impact
            score = max(raw_score, min_score)

        scores[mutation] = round(score, 3)
    
    return scores