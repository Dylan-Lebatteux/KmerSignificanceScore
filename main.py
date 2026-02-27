#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
K-mer Significance Score (KSS) Analysis Tool

This script performs viral genomic sequence analysis using k-mer based scoring
that integrates discriminative, mutational, and protein-level features.

Usage:
    python main.py [config_files...]

    Examples:
        python main.py data/Human_betaherpesvirus_5/config.yaml
        python main.py data/*/config.yaml
        python main.py  # Uses all config.yaml files in data/

Configuration is loaded from YAML files (one per virus dataset).
"""

import os
import sys
import time
import io

# Fix Windows console encoding
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Import KSS modules
# Note: sys.path modification allows running main.py from project root
sys.path.insert(0, 'src')
from src import utils
from src import kss
from src import kanalyzer
from src import pipeline


# ============================================================================
# CONFIGURATION LOADING (now handled by pipeline module)
# ============================================================================

# Note: Configuration loading functions have been moved to src/pipeline.py
# This eliminates code duplication between main.py and recalculate_scores.py


# ============================================================================
# MAIN PROCESSING FUNCTIONS
# ============================================================================

def process_single_gene(gene: str, dataset_name: str, info: dict, parameters: dict,
                       gene_num: int, total_genes: int) -> None:
    """
    Process a single gene through the complete KSS pipeline.

    This function processes one gene at a time to optimize memory usage,
    especially important for large datasets (e.g., HIV-1 with 12,000+ sequences).

    Args:
        gene: Gene name to process
        dataset_name: Name of the dataset being processed
        info: Dictionary containing input/output folder paths
        parameters: Dictionary of analysis parameters
        gene_num: Current gene number (for progress display)
        total_genes: Total number of genes to process

    Returns:
        None. Results are saved to gene-specific folder.
    """
    print(f"\n[Gene {gene_num}/{total_genes}] Processing {gene}...")
    gene_start_time = time.time()

    # Create gene-specific info (only analyze this gene)
    gene_info = info.copy()
    gene_info['cds_selection'] = gene

    try:
        # Step 1: Analyze sequences and extract k-mers for this gene only
        print(f"  [1/3] Analyzing {gene} sequences and extracting k-mers...")
        results = kanalyzer.analyze_records(gene_info, parameters)

        # Step 2: Compile results and calculate KSS scores for this gene
        print(f"  [2/3] Compiling and calculating scores for {gene}...")
        # OPTIMIZATION: Pass target_gene to process only this gene (instead of iterating over all genes)
        compiled_results = kss.compile_results(results, parameters, target_gene=gene)

        compiled_results, discriminative_scores = kss.compute_kss_scores(
            compiled_results,
            results,
            gene_info,
            parameters,
            target_gene=gene  # OPTIMIZATION: Pass target_gene to process only this gene
        )

        # Step 3: Save results for this gene
        print(f"  [3/3] Saving {gene} results...")
        gene_output_dir = os.path.join(info['output_folder'], gene, "results")
        os.makedirs(gene_output_dir, exist_ok=True)

        saved_files = 0

        # Save raw results if requested
        if parameters.get("save_raw", False):
            gene_results = {
                "_metadata": results.get("_metadata", {}),
                "genes": {gene: results["genes"][gene]} if gene in results.get("genes", {}) else {}
            }
            raw_path = os.path.join(gene_output_dir, f"{gene}_results.json")
            utils.save_data_as_json(gene_results, raw_path)
            saved_files += 1

        # Save compiled results (includes all scores: protein, mutational, discriminative, KSS)
        gene_compiled = {gene: compiled_results[gene]} if gene in compiled_results else {}
        compiled_path = os.path.join(gene_output_dir, f"{gene}_compiled_results.json")
        utils.save_data_as_json(gene_compiled, compiled_path)
        saved_files += 1

        # Note: discriminative_scores are already included in compiled_results
        # No need to save them separately

        gene_time = time.time() - gene_start_time
        print(f"  ✓ {gene}: {saved_files} file(s) saved | Time: {gene_time:.1f}s")

    except Exception as e:
        import traceback
        print(f"  ✗ Error processing {gene}:")
        print(f"  {str(e)}")
        traceback.print_exc()
        raise


def process_dataset(dataset_name: str, info: dict, parameters: dict) -> None:
    """
    Process a specific viral dataset using K-mer Significance Score analysis.

    This function processes genes sequentially (gene-by-gene) to minimize
    memory usage and provide better progress tracking. Each gene is fully
    processed (analyzed, scored, saved) before moving to the next.

    Args:
        dataset_name: Name of the dataset being processed
        info: Dictionary containing:
            - input_folder: Path to input data directory
            - output_folder: Base path for output (gene-specific subdirs created)
            - cds_selection: Comma-separated list of genes/CDS to analyze
        parameters: Dictionary of analysis parameters

    Returns:
        None. Results are saved to gene-specific folders.
    """
    print(f"\n{'='*70}")
    print(f"Processing dataset: {dataset_name}")
    print(f"{'='*70}")

    genes = info['cds_selection'].split(',')
    total_genes = len(genes)
    print(f"Genes to process: {total_genes} ({info['cds_selection']})")
    print(f"Strategy: Sequential gene-by-gene processing (optimized for memory)")

    try:
        dataset_start_time = time.time()

        # Process each gene sequentially
        for idx, gene in enumerate(genes, 1):
            process_single_gene(gene, dataset_name, info, parameters, idx, total_genes)

        total_time = time.time() - dataset_start_time

        print(f"\n{'='*70}")
        print(f"✓ Completed: {dataset_name}")
        print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
        print(f"Average time per gene: {total_time/total_genes:.1f}s")
        print(f"{'='*70}\n")

    except Exception as e:
        import traceback
        print(f"\n{'!'*70}")
        print(f"✗ Error processing {dataset_name}:")
        print(f"  {str(e)}")
        print("\nFull traceback:")
        traceback.print_exc()
        print(f"{'!'*70}\n")
        raise


def main():
    """
    Main execution function for K-mer Significance Score analysis.

    Loads configuration from YAML files and processes viral datasets.
    """
    print("\n" + "="*70)
    print("K-mer Significance Score (KSS) Analysis")
    print("Viral Genomic Sequence Analysis Tool")
    print("="*70)

    # Find and load configuration files
    config_files = pipeline.find_config_files(sys.argv[1:])

    print(f"\nFound {len(config_files)} configuration file(s):")
    for config_file in config_files:
        print(f"  • {config_file}")

    # Load all configurations
    datasets = []
    for config_file in config_files:
        try:
            dataset_name, dataset_info, parameters = pipeline.load_config(config_file)
            datasets.append((dataset_name, dataset_info, parameters, config_file))
            genes = dataset_info['cds_selection'].split(',')
            print(f"    ✓ {dataset_name}: {len(genes)} genes, k={parameters['k']}")
        except Exception as e:
            print(f"    ✗ Error loading {config_file}: {e}")
            continue

    if not datasets:
        print("\nError: No valid configurations loaded!")
        sys.exit(1)

    # Process each dataset
    successful = 0
    failed = 0

    for idx, (dataset_name, info, parameters, config_file) in enumerate(datasets, 1):
        print(f"\n[Dataset {idx}/{len(datasets)}] {dataset_name}")
        try:
            process_dataset(dataset_name, info, parameters)
            successful += 1
        except Exception:
            failed += 1
            continue

    # Summary
    print("\n" + "="*70)
    print("Analysis Summary")
    print("="*70)
    print(f"  Successful: {successful}/{len(datasets)}")
    if failed > 0:
        print(f"  Failed: {failed}/{len(datasets)}")
    print("="*70)

    if failed == 0:
        print("\n✓ All analyses completed successfully!\n")
    else:
        print(f"\n✗ {failed} dataset(s) failed. Check error messages above.\n")


if __name__ == "__main__":
    main()
