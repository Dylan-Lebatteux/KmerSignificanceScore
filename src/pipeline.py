#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
KSS Pipeline Module

This module provides reusable workflow functions for the K-mer Significance Score
analysis pipeline. It centralizes common operations such as configuration loading,
score recalculation, and dataset processing to promote code reuse and maintainability.

Main functions:
    find_config_files: Find configuration files from arguments or defaults
    load_config: Load and parse configuration from YAML files
    recalculate_scores_for_gene: Recalculate scores for a single gene
    recalculate_scores_from_config: Recalculate scores for entire dataset
"""

import os
import sys
import glob
import yaml
from typing import Dict, Any, Tuple, Optional, List

# Import KSS modules
try:
    from . import utils
    from . import kss
except ImportError:
    # Fallback for running as standalone script
    import utils
    import kss


def find_config_files(args: List[str]) -> List[str]:
    """
    Find configuration files from command line arguments or default locations.

    Args:
        args: Command line arguments (e.g., sys.argv[1:])

    Returns:
        List of configuration file paths

    Raises:
        SystemExit: If no configuration files are found

    Notes:
        - Supports glob patterns (e.g., "data/*/config.yaml")
        - If no arguments provided, searches for data/*/config.yaml
        - Prints usage information and exits if no files found

    Examples:
        >>> files = find_config_files(['data/CMV/config.yaml'])
        >>> files = find_config_files(['data/*/config.yaml'])
        >>> files = find_config_files([])  # Uses default pattern
    """
    config_files = []

    if args:
        # Use provided config files or patterns
        for arg in args:
            matched = glob.glob(arg)
            config_files.extend(matched)
    else:
        # Default: find all config.yaml in data/ subdirectories
        config_files = glob.glob("data/*/config.yaml")

    if not config_files:
        print("Error: No configuration files found!")
        print("\nUsage:")
        print("  python main.py [config_files...]")
        print("\nExamples:")
        print("  python main.py data/Human_betaherpesvirus_5/config.yaml")
        print("  python main.py data/*/config.yaml")
        print("  python main.py  # Uses all config.yaml files in data/")
        sys.exit(1)

    return config_files


def load_config(config_path: str) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """
    Load configuration from YAML file.

    Args:
        config_path: Path to YAML configuration file

    Returns:
        Tuple of (dataset_name, dataset_info, parameters)
            - dataset_name: Name of the dataset
            - dataset_info: Dictionary with input_folder, output_folder, cds_selection
            - parameters: Dictionary with all analysis parameters

    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config file is malformed
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Configuration file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Extract dataset info
    dataset = config['dataset']
    dataset_name = dataset['name']

    # Build dataset info dictionary
    base_folder = os.path.dirname(config_path)
    dataset_info = {
        "input_folder": base_folder,
        "output_folder": base_folder,
        "cds_selection": ",".join(dataset['genes'])
    }

    # Extract parameters
    params = config['parameters']
    parameters = {
        "k": params['k'],
        "threshold": params['threshold'],
        "discriminative_weight": params['weights']['discriminative'],
        "mutational_weight": params['weights']['mutational'],
        "protein_weight": params['weights']['protein'],
        "mutational_matrix": params['scoring']['mutational_matrix'],
        "substitution_matrix": params['scoring']['substitution_matrix'],
        "indel_score": params['scoring'].get('indel_score', 1.0),
        "open_gap_score": params['alignment']['open_gap_score'],
        "extend_gap_score": params['alignment']['extend_gap_score'],
        "save_raw": params['output'].get('save_raw', False)
    }

    return dataset_name, dataset_info, parameters


def recalculate_scores_for_gene(gene: str,
                                base_folder: str,
                                parameters: Dict[str, Any],
                                verbose: bool = True) -> bool:
    """
    Recalculate KSS scores for a single gene from existing raw results.

    This function loads pre-computed raw results (mutations and alignments)
    and recalculates only the compiled results and KSS scores. It does not
    re-run sequence alignments or mutation detection.

    Args:
        gene: Gene name to process
        base_folder: Base dataset folder path
        parameters: Analysis parameters dictionary
        verbose: If True, prints progress information

    Returns:
        True if successful, False if raw results not found or error occurred

    Notes:
        - Requires existing {gene}_results.json file
        - Overwrites existing {gene}_compiled_results.json
        - Uses same parameters as original analysis
        - Much faster than full analysis (no alignments)
    """
    if verbose:
        print(f"\nProcessing gene: {gene}")

    # Load existing raw results
    results_path = os.path.join(base_folder, gene, "results", f"{gene}_results.json")

    if not os.path.exists(results_path):
        if verbose:
            print(f"  ✗ Raw results not found: {results_path}")
        return False

    if verbose:
        print(f"  Loading raw results: {results_path}")

    try:
        results = utils.load_data_from_json(results_path)
    except Exception as e:
        if verbose:
            print(f"  ✗ Error loading raw results: {e}")
        return False

    # Prepare info dict
    infos = {
        "input_folder": base_folder,
        "output_folder": base_folder
    }

    try:
        # Step 1: Compile results
        if verbose:
            print(f"  [1/2] Compiling results...")
        compiled_results = kss.compile_results(results, parameters, target_gene=gene)

        # Step 2: Compute scores
        if verbose:
            print(f"  [2/2] Computing KSS scores...")
        compiled_results, discriminative_scores = kss.compute_kss_scores(
            compiled_results,
            results,
            infos,
            parameters,
            target_gene=gene
        )

        # Save compiled results
        gene_output_dir = os.path.join(base_folder, gene, "results")
        os.makedirs(gene_output_dir, exist_ok=True)

        compiled_path = os.path.join(gene_output_dir, f"{gene}_compiled_results.json")

        gene_compiled = {gene: compiled_results[gene]} if gene in compiled_results else {}
        utils.save_data_as_json(gene_compiled, compiled_path)

        if verbose:
            print(f"  ✓ Saved: {compiled_path}")

        return True

    except Exception as e:
        if verbose:
            import traceback
            print(f"  ✗ Error recalculating scores: {e}")
            traceback.print_exc()
        return False


def recalculate_scores_from_config(config_path: str,
                                   target_gene: Optional[str] = None,
                                   verbose: bool = True) -> Dict[str, bool]:
    """
    Recalculate KSS scores for dataset from configuration file.

    Main entry point for score recalculation workflow. Loads configuration
    and recalculates scores for all or specified genes.

    Args:
        config_path: Path to YAML configuration file
        target_gene: Optional specific gene to process. If None, processes all genes
        verbose: If True, prints detailed progress information

    Returns:
        Dictionary mapping gene names to success status (True/False)

    Example:
        >>> results = recalculate_scores_from_config(
        ...     "data/Human_betaherpesvirus_5/config.yaml",
        ...     target_gene="UL55"
        ... )
        >>> print(results)
        {'UL55': True}

    Notes:
        - Only recalculates genes with existing raw results
        - Returns partial success if some genes fail
        - Existing compiled results will be overwritten
    """
    # Load configuration
    if verbose:
        print(f"\nLoading configuration: {config_path}")

    try:
        dataset_name, dataset_info, parameters = load_config(config_path)
    except Exception as e:
        if verbose:
            print(f"✗ Error loading configuration: {e}")
        return {}

    base_folder = dataset_info['input_folder']

    # Determine genes to process
    all_genes = dataset_info['cds_selection'].split(',')
    genes_to_process = [target_gene] if target_gene else all_genes

    if verbose:
        print(f"\n{'='*70}")
        print(f"Recalculating KSS scores: {dataset_name}")
        print(f"{'='*70}")
        print(f"Genes to process: {', '.join(genes_to_process)}")
        print(f"{'='*70}")

    # Recalculate scores for each gene
    results = {}
    successful = 0
    failed = 0

    for gene in genes_to_process:
        success = recalculate_scores_for_gene(gene, base_folder, parameters, verbose=verbose)
        results[gene] = success

        if success:
            successful += 1
        else:
            failed += 1

    # Summary
    if verbose:
        print(f"\n{'='*70}")
        print(f"Summary")
        print(f"{'='*70}")
        print(f"  Successful: {successful}/{len(genes_to_process)}")
        if failed > 0:
            print(f"  Failed: {failed}/{len(genes_to_process)}")
        print(f"{'='*70}\n")

    return results
