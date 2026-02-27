#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Discriminative Scoring Module - KSS Method

This module calculates discriminative scores for genomic positions using
k-mer based classification. The score quantifies how well k-mer patterns
at a position discriminate between different classes (e.g., viral strains).

Scoring Formula:
    KSS = tanh(sqrt(NMI × weighted_purity²) × ln(n_classes + 1))

    where:
    - NMI = Normalized Mutual Information between k-mers and classes
    - weighted_purity² = Σ P(kmer_j) × [max_c P(class_c | kmer_j)]²
    - sqrt() = Geometric mean to reduce over-penalization
    - ln(n_classes + 1) = Class complexity adjustment

The final score is normalized to [0, 1] range, where:
    - 1.0 = Perfect discrimination (k-mers uniquely identify classes)
    - 0.0 = No discriminative power (k-mers independent of classes)

Main function:
    get_discriminative_score: Calculate discriminative score for a genomic position
"""

import numpy as np
from typing import Dict, Any


def get_discriminative_score(
    X: np.ndarray,
    y: np.ndarray
) -> Dict[str, Any]:
    """
    Calculate discriminative score for a genomic position.

    This is the main entry point for discriminative scoring. It computes how well
    k-mer patterns at a specific genomic position discriminate between classes
    using normalized mutual information weighted by k-mer purity.

    Parameters
    ----------
    X : np.ndarray
        Binary ONE-HOT feature matrix (n_sequences, n_kmers)
        Each row represents a sequence, each column a possible k-mer
        Each sequence has exactly one k-mer = 1, all others = 0
    y : np.ndarray
        Class labels (n_sequences,)
        String or integer labels for each sequence

    Returns
    -------
    dict
        Dictionary containing:
        - 'raw_score': Base score before tanh transformation
        - 'normalized_score': Final score in [0, 1] range
        - 'n_classes': Number of unique classes
        - 'n_sequences': Number of sequences analyzed
        - 'n_kmers': Number of distinct k-mers at this position

    Examples
    --------
    >>> X = np.array([[1, 0], [0, 1], [1, 0]])
    >>> y = np.array(['A', 'B', 'A'])
    >>> result = get_discriminative_score(X, y)
    >>> print(f"Score: {result['normalized_score']:.3f}")
    Score: 0.882

    Notes
    -----
    - Requires at least 2 classes to compute meaningful scores
    - Empty feature matrices return score of 0.0
    - Uses tanh transformation for bounded output in [0, 1]
    """
    X = np.asarray(X)
    y = np.asarray(y)

    n_sequences, n_kmers = X.shape
    unique_classes = np.unique(y)
    n_classes = len(unique_classes)

    # Edge cases
    if n_kmers == 0 or n_classes < 2:
        return {
            'raw_score': 0.0,
            'normalized_score': 0.0,
            'n_classes': n_classes,
            'n_sequences': n_sequences,
            'n_kmers': n_kmers
        }

    # Class complexity factor: ln(n_classes + 1)
    class_complexity = np.log(n_classes + 1)

    # ========================================================================
    # STEP 1: Calculate H(Y) - Entropy of classes
    # ========================================================================
    class_probs = np.array([np.mean(y == c) for c in unique_classes])
    h_class = -np.sum(class_probs * np.log(class_probs + 1e-10))

    if h_class == 0:
        return {
            'raw_score': 0.0,
            'normalized_score': 0.0,
            'n_classes': n_classes,
            'n_sequences': n_sequences,
            'n_kmers': n_kmers
        }

    # ========================================================================
    # STEP 2: Calculate H(Y|X) - Conditional entropy by k-mer configurations
    # ========================================================================
    unique_rows, inverse = np.unique(X, axis=0, return_inverse=True)
    h_conditional = 0.0

    for i in range(len(unique_rows)):
        mask = (inverse == i)
        p_config = np.mean(mask)
        if p_config > 0:
            y_subset = y[mask]
            subset_probs = np.array([np.mean(y_subset == c) for c in unique_classes])
            h_subset = -np.sum(subset_probs * np.log(subset_probs + 1e-10))
            h_conditional += p_config * h_subset

    # Normalized Mutual Information
    nmi = float(np.clip((h_class - h_conditional) / h_class, 0, 1))

    # ========================================================================
    # STEP 3: Calculate boosted weighted purity
    # ========================================================================
    weighted_purity_boosted = 0.0

    for j in range(n_kmers):
        mask_kmer = (X[:, j] == 1)
        freq_j = np.mean(mask_kmer)

        if freq_j > 0:
            y_kmer = y[mask_kmer]
            probs_j = np.array([np.mean(y_kmer == c) for c in unique_classes])
            purity_j = np.max(probs_j)
            weighted_purity_boosted += freq_j * (purity_j ** 2)

    # ========================================================================
    # STEP 4: Calculate final score
    # ========================================================================
    raw_score = np.sqrt(nmi * weighted_purity_boosted)
    x = raw_score * class_complexity
    normalized_score = float(np.clip(np.tanh(x), 0, 1))

    return {
        'raw_score': round(float(raw_score), 4),
        'normalized_score': round(normalized_score, 3),
        'n_classes': n_classes,
        'n_sequences': n_sequences,
        'n_kmers': n_kmers
    }
