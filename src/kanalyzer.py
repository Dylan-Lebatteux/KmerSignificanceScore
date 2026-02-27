#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Imports
import os
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import BiopythonWarning
from Bio.Align import PairwiseAligner, substitution_matrices
from . import mutation_score

# Suppress BiopythonWarning
warnings.filterwarnings("ignore", category=BiopythonWarning)


def initialize_aligner(parameters: dict) -> PairwiseAligner:
    """
    Initialize a pairwise sequence aligner with specified parameters.

    Args:
        parameters: Dictionary containing alignment parameters:
            - substitution_matrix: Name of substitution matrix (e.g., 'BLOSUM62')
            - open_gap_score: Penalty for opening a gap
            - extend_gap_score: Penalty for extending a gap

    Returns:
        Configured PairwiseAligner object for global sequence alignment
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load(parameters["substitution_matrix"])
    aligner.open_gap_score = parameters["open_gap_score"]
    aligner.extend_gap_score = parameters["extend_gap_score"]

    return aligner


def align_nucleotides(protein_a: str, protein_b: str, nucleotide_a: str, nucleotide_b: str) -> tuple:
    """
    Align nucleotide sequences based on aligned protein sequences.

    This function takes protein alignments and back-translates them to
    nucleotide alignments, preserving codon structure (3 nucleotides per amino acid).

    Args:
        protein_a: First aligned protein sequence (may contain gaps '-')
        protein_b: Second aligned protein sequence (may contain gaps '-')
        nucleotide_a: Nucleotide sequence corresponding to protein_a
        nucleotide_b: Nucleotide sequence corresponding to protein_b

    Returns:
        Tuple of (aligned_nucleotide_a, aligned_nucleotide_b) as strings
    """
    nucleotide_a = Seq(nucleotide_a)
    nucleotide_b = Seq(nucleotide_b)

    aligned_nuc_a = []
    aligned_nuc_b = []
    index_a = index_b = 0

    for aa_a, aa_b in zip(protein_a, protein_b):
        if aa_a == '-':
            aligned_nuc_a.append('---')
        else:
            aligned_nuc_a.append(str(nucleotide_a[index_a:index_a+3]))
            index_a += 3

        if aa_b == '-':
            aligned_nuc_b.append('---')
        else:
            aligned_nuc_b.append(str(nucleotide_b[index_b:index_b+3]))
            index_b += 3

    return ''.join(aligned_nuc_a), ''.join(aligned_nuc_b)


def replace_first_three_gaps(sequence: str, replacement: str = '') -> str:
    """
    Replace the first three gap characters '-' in a sequence with a replacement string.

    Used for handling insertions in nucleotide alignment (one codon = 3 nucleotides).

    Args:
        sequence: Input sequence containing gap characters
        replacement: String to replace gaps with (default: empty string)

    Returns:
        Modified sequence with first three gaps replaced
    """
    count = 0
    modified_sequence = []
    for char in sequence:
        if char == '-' and count < 3:
            modified_sequence.append(replacement)
            count += 1
        else:
            modified_sequence.append(char)
    return ''.join(modified_sequence)


def identify_mutations(infos: dict, parameters: dict, ref_nuc_seq: str,
                      ref_aa_seq: str, query_nuc_seq: str, query_aa_seq: str, aligner=None) -> list:
    """
    Identify nucleotide and amino acid mutations between reference and query sequences.

    This function performs k-mer based mutation analysis by:
    1. Aligning protein sequences
    2. Back-translating to nucleotide alignments
    3. Extracting k-mers and identifying variations
    4. Annotating amino acid changes

    Only positions with actual mutations are returned (incremental storage).

    Args:
        infos: Dataset information dictionary
        parameters: Analysis parameters including k-mer size ('k')
        ref_nuc_seq: Reference nucleotide sequence
        ref_aa_seq: Reference amino acid sequence
        query_nuc_seq: Query nucleotide sequence to compare
        query_aa_seq: Query amino acid sequence to compare

    Returns:
        List of mutation objects, each containing:
        - position: Nucleotide position (1-indexed)
        - ref_kmer: Reference k-mer sequence
        - alt_kmer: Alternative k-mer sequence
        - aa_changes: List of amino acid changes with impact scores
    """
    mutations = []
    k = parameters["k"]

    # Perform sequence alignments
    if aligner is None:
        aligner = initialize_aligner(parameters)
    alignments = aligner.align(ref_aa_seq, query_aa_seq)
    aligned_ref_aa, aligned_query_aa = alignments[0]
    aligned_ref_nuc, aligned_query_nuc = align_nucleotides(aligned_ref_aa, aligned_query_aa, ref_nuc_seq, query_nuc_seq)

    n_insertion = 0  # Track insertions

    # Analyze mutations - only store if different from reference
    for i in range(0, len(aligned_ref_nuc) - k + 1, k):
        adjusted_k = k
        current_pos = i + (n_insertion * 3)

        # Get current k-mers
        ref_kmer = aligned_ref_nuc[current_pos:current_pos + adjusted_k]
        query_kmer = aligned_query_nuc[current_pos:current_pos + adjusted_k]

        # Handle gaps in reference sequence
        while '-' in ref_kmer and current_pos + adjusted_k < len(aligned_ref_nuc):
            next_block = current_pos + adjusted_k
            ref_kmer = replace_first_three_gaps(ref_kmer) + aligned_ref_nuc[next_block:next_block + 3]
            adjusted_k += 3
            query_kmer = aligned_query_nuc[current_pos:current_pos + adjusted_k]

        # Update insertion tracking
        current_insertion = (adjusted_k - k) // 3
        n_insertion += current_insertion

        # Skip if k-mers are identical (incremental storage)
        if ref_kmer == query_kmer:
            continue

        # Calculate amino acid positions
        aa_start = current_pos // 3
        aa_end = min(aa_start + (adjusted_k // 3), len(aligned_ref_aa))

        ref_aa_kmer = aligned_ref_aa[aa_start:aa_end]
        query_aa_kmer = aligned_query_aa[aa_start:aa_end]

        # Process amino acid changes
        aa_changes = []
        aa_position = aa_start - n_insertion
        processed_insertions = 0

        for idx, (ref_aa, query_aa) in enumerate(zip(ref_aa_kmer, query_aa_kmer)):
            if ref_aa != query_aa:
                # Calculate mutation position
                mut_pos = max(0, aa_position + idx + 1 + current_insertion - processed_insertions)

                # Create mutation annotation
                notation = f"{ref_aa}{mut_pos}{query_aa}"
                if ref_aa == '-':
                    processed_insertions += 1

                # Handle duplicate positions with suffix
                original_notation = notation
                suffix = ""
                while any(change["notation"] == notation for change in aa_changes):
                    suffix += "*"
                    notation = original_notation + suffix

                aa_changes.append({
                    "notation": notation,
                    "mut_score": 0  # Placeholder, will be updated below
                })

        # Store mutation only if there are changes
        if aa_changes or ref_kmer != query_kmer:
            mutations.append({
                "position": i + 1,
                "ref_kmer": ref_kmer,
                "alt_kmer": query_kmer,
                "aa_changes": aa_changes
            })

    # Calculate mutation scores for all amino acid changes using substitution matrix
    if mutations:
        # Collect all unique mutation notations
        all_notations = {}
        for mutation in mutations:
            for aa_change in mutation["aa_changes"]:
                notation = aa_change["notation"]
                # Remove suffix for score calculation (e.g., "A10G*" -> "A10G")
                clean_notation = notation.rstrip("*")
                if clean_notation not in all_notations:
                    all_notations[clean_notation] = 0

        # Get mutation scores from substitution matrix
        mutation_scores = mutation_score.get_mutational_scores(
            all_notations, parameters["mutational_matrix"]
        )

        # Update mutation scores in aa_changes
        for mutation in mutations:
            for aa_change in mutation["aa_changes"]:
                notation = aa_change["notation"]
                clean_notation = notation.rstrip("*")
                if clean_notation in mutation_scores:
                    aa_change["mut_score"] = mutation_scores[clean_notation]

    return mutations


def calculate_similarity(seq_a: str, seq_b: str) -> float:
    """
    Calculate the percentage similarity between two sequences.

    Args:
        seq_a: First sequence
        seq_b: Second sequence

    Returns:
        Percentage similarity (0-100)
    """
    matches = sum(a == b for a, b in zip(seq_a, seq_b))
    total = len(seq_a)
    return (100 * matches / total) if total else 0


def analyze_records(infos: dict, parameters: dict) -> dict:
    """
    Analyze genomic records from GenBank and FASTA files to identify mutations.

    This is the main entry point for sequence analysis. It:
    1. Loads reference sequences from GenBank files
    2. Loads query sequences from FASTA files
    3. Identifies mutations for each query sequence

    Uses incremental storage: only positions with mutations are stored.

    Args:
        infos: Dictionary containing:
            - input_folder: Path to data directory
            - cds_selection: Comma-separated list of genes/CDS to analyze
        parameters: Analysis parameters (k-mer size, alignment parameters, etc.)

    Returns:
        Dictionary structure:
        {
            "_metadata": {version, k, format, ...},
            "genes": {
                gene: {
                    "metadata": {reference_length, genbank_id, ...},
                    "sequences": {
                        sequence_id: {
                            "class": class_label,
                            "mutations": [...]
                        }
                    }
                }
            }
        }
    """
    # Initialize results structure with metadata
    results = {
        "_metadata": {
            "version": "2.0.0",
            "tool": "K-mer Significance Score (KSS)",
            "k": parameters["k"]
        },
        "genes": {}
    }

    cds_list = infos["cds_selection"].split(",")
    total_sequences = 0

    for cds in cds_list:
        # Get reference data
        gb_dir = os.path.join(infos["input_folder"], cds)
        gb_files = [f for f in os.listdir(gb_dir) if f.lower().endswith('.gb')]
        if not gb_files:
            print(f"  ⚠ No GenBank file found for {cds}")
            continue

        gb_path = os.path.join(gb_dir, gb_files[0])
        gb_record = list(SeqIO.parse(gb_path, "genbank"))[0]

        # Extract reference sequences
        ref_nuc = ""
        ref_aa = ""
        for feature in gb_record.features:
            if feature.type == "CDS" and feature.qualifiers.get("gene", [""])[0] == cds:
                ref_nuc = str(feature.location.extract(gb_record.seq))
                ref_aa = feature.qualifiers.get("translation", [""])[0]
                break

        if not ref_nuc or not ref_aa:
            print(f"  ⚠ CDS {cds} not found in GenBank file")
            continue

        # Get query data
        fasta_dir = os.path.join(infos["input_folder"], cds)
        fasta_files = [f for f in os.listdir(fasta_dir) if f.lower().endswith('.fasta')]
        if not fasta_files:
            print(f"  ⚠ No FASTA file found for {cds}")
            continue

        fasta_path = os.path.join(fasta_dir, fasta_files[0])

        # Count sequences first without loading them all
        print(f"  Counting {cds} sequences...", end='', flush=True)
        num_sequences = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
        print(f" {num_sequences} sequences found")

        # Initialize gene structure
        k = parameters["k"]
        num_positions = (len(ref_nuc) - k + 1 + k - 1) // k

        results["genes"][cds] = {
            "metadata": {
                "reference_length": len(ref_nuc),
                "reference_aa_length": len(ref_aa),
                "genbank_id": gb_record.id,
                "total_positions": num_positions,
                "num_sequences": num_sequences
            },
            "sequences": {}
        }

        # Pre-initialize aligner once (reuse across all sequences)
        aligner = initialize_aligner(parameters)

        # Process each query sequence - stream instead of loading all at once
        print(f"  Processing {cds} sequences: ", end='', flush=True)

        for idx, record in enumerate(SeqIO.parse(fasta_path, "fasta"), 1):
            # Progress indicator every 100 sequences
            if idx % 100 == 0 or idx == num_sequences:
                percent = (idx / num_sequences) * 100
                print(f"\r  Processing {cds} sequences: {idx}/{num_sequences} ({percent:.1f}%)", end='', flush=True)
            query_id = record.id
            query_nuc = str(record.seq)
            # Translate nucleotide sequence to amino acids
            # Note: Biopython's PairwiseAligner does not support 'J' (ambiguous Leu/Ile)
            # We normalize 'J' to 'L' (Leucine) as the closest standard amino acid
            # This occurs rarely and only affects sequences with non-standard translation
            query_aa = str(Seq(query_nuc).translate()).replace('J', 'L')

            # Extract class label from sequence header (format: >ID|Virus|Gene|CLASS)
            class_label = "unknown"
            if '|' in query_id:
                parts = query_id.split('|')
                # Class is the last field
                if len(parts) > 0:
                    class_label = parts[-1]

            # Identify mutations (returns list of mutations only)
            mutations = identify_mutations(
                infos, parameters, ref_nuc, ref_aa, query_nuc, query_aa, aligner
            )

            results["genes"][cds]["sequences"][query_id] = {
                "class": class_label,
                "mutations": mutations
            }

            total_sequences += 1

        # Print summary for this gene
        print(f"\n  ✓ {cds}: {num_sequences} sequences processed")

    # Update global metadata
    results["_metadata"]["sequences_analyzed"] = total_sequences

    return results