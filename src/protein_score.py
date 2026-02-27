#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Protein Functional Importance Scoring Module

This module calculates protein functional importance scores based on UniProt
annotations and literature evidence. It integrates multiple biological evidence
types to quantify how well-characterized and functionally important a protein is.

Scoring Components (total: 100 points):
    - Protein existence evidence: 15 points
    - Molecular function annotations: 20 points
    - Biological process annotations: 20 points
    - Cellular component annotations: 5 points
    - Structural annotations: 10 points
    - Protein-protein interactions: 5 points
    - Post-translational modifications: 5 points
    - 3D structure availability: 5 points
    - Drug target information: 2.5 points
    - Publication references: 12.5 points

The final score is normalized to [0, 1] range, where:
    - 1.0 = Extensively characterized protein with strong evidence
    - 0.0 = Poorly characterized or hypothetical protein

Main functions:
    get_protein_score: Calculate overall protein importance score
    fetch_uniprot_data: Retrieve and cache UniProt entries
    calculate_protein_score: Compute score from cached UniProt JSON
    score_* functions: Individual scoring components
"""

import os
import json
import time
from io import StringIO
from typing import Dict, Any, Optional, Callable
import requests
import pandas as pd


def _retry_request(request_func: Callable, max_retries: int = 3,
                   initial_delay: float = 1.0, verbose: bool = False) -> requests.Response:
    """
    Retry an HTTP request with exponential backoff.

    Args:
        request_func: Callable that returns a requests.Response
        max_retries: Maximum number of retry attempts (default: 3)
        initial_delay: Initial delay in seconds between retries (default: 1.0)
        verbose: If True, print retry information

    Returns:
        requests.Response object if successful

    Raises:
        requests.exceptions.RequestException: If all retries fail
    """
    retry_delay = initial_delay

    for attempt in range(max_retries):
        try:
            response = request_func()
            response.raise_for_status()
            return response
        except (requests.exceptions.RequestException, requests.exceptions.Timeout) as e:
            if attempt < max_retries - 1:
                if verbose:
                    print(f"Request failed (attempt {attempt + 1}/{max_retries}), "
                          f"retrying in {retry_delay}s...")
                time.sleep(retry_delay)
                retry_delay *= 2  # Exponential backoff
            else:
                raise  # Re-raise on last attempt


def fetch_uniprot_data(taxon_id: str,
                      gene_name: str,
                      verbose: bool = False) -> str:
    """
    Retrieve UniProt data for a given taxon ID and gene name.

    This function queries the UniProt REST API to find reviewed entries matching
    the specified taxonomy and gene name. It selects the best-annotated entry
    based on number of publications and caches the JSON data locally.

    Args:
        taxon_id: NCBI Taxonomy ID of the organism (e.g., "9606" for human)
        gene_name: Name of the gene to query. Must match exactly in UniProt's
                  Gene Names field. Special cases:
                  - "ORF1ab" is converted to "rep"
                  - "pol" is converted to "gag-pol"
        verbose: If True, prints detailed information about retrieval process

    Returns:
        UniProt entry ID (accession) if successful, or an error/status message
        string starting with "No" or "Error" if unsuccessful

    Notes:
        - Only searches reviewed (Swiss-Prot) entries
        - Ranks results by number of PubMed publications
        - Downloads and caches full JSON data in 'uniprot/' subdirectory
        - Cache prevents redundant API calls
        - Returns first (best) matching entry
    """
    try:
        # Use requests directly with UniProt REST API
        if gene_name == "ORF1ab":
            gene_name = "rep"
        elif gene_name == "pol":
            gene_name = "gag-pol"

        # Query UniProt API with retry logic
        response = _retry_request(
            lambda: requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query": f"{taxon_id} {gene_name} AND reviewed:true",
                    "format": "tsv",
                    "fields": "accession,id,gene_names,annotation_score,lit_pubmed_id",
                    "size": 100
                },
                timeout=30
            ),
            verbose=verbose
        )
        
        # Convert response to DataFrame
        df = pd.read_csv(StringIO(response.text), sep='\t')
        
        if verbose:
            print(f"Retrieved {len(df)} entries for taxon {taxon_id}")
        
        # Filter to keep only entries with exact gene name match
        if 'Gene Names' in df.columns:
            # Create explicit copy to avoid SettingWithCopyWarning
            df_filtered = df[df['Gene Names'].apply(
                lambda x: gene_name in str(x).split() if pd.notna(x) else False
            )].copy()
            
            if verbose:
                print(f"Filtered to {len(df_filtered)} entries matching gene name '{gene_name}'")
        else:
            df_filtered = df.copy()
            if verbose:
                print("Warning: 'Gene Names' column not found")
        
        if df_filtered.empty:
            return f"No entries found with exact gene name '{gene_name}' in taxon {taxon_id}"
        
        # Count and sort by number of PubMed publications
        if 'PubMed ID' in df_filtered.columns:
            # Count number of publications
            df_filtered['PubMed Count'] = df_filtered['PubMed ID'].apply(
                lambda x: len(str(x).split(';')) if pd.notna(x) and str(x).strip() else 0
            )
            # Sort by number of publications (descending)
            df_filtered = df_filtered.sort_values(by=['PubMed Count'], ascending=False)
            # Remove PubMed ID column after calculation
            df_filtered = df_filtered.drop(columns=['PubMed ID'])
            
            if verbose:
                print("Sorted entries by number of PubMed publications")
        
        # Display filtered and sorted results
        if verbose:
            print(df_filtered)
        
        if df_filtered.empty:
            return f"No suitable entries found after filtering for gene {gene_name}"
            
        # Select entry with most publications
        entry_col = 'Entry' if 'Entry' in df_filtered.columns else 'Accession'
        entry_id = df_filtered.iloc[0][entry_col]
        
        # Cache handling
        cache_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'uniprot')
        os.makedirs(cache_dir, exist_ok=True)
        json_file_path = os.path.join(cache_dir, f'{entry_id}.json')
        
        # Download and cache if not exists
        if not os.path.exists(json_file_path):
            if verbose:
                print(f"Downloading JSON for {entry_id}...")

            # Download JSON with retry logic
            json_response = _retry_request(
                lambda: requests.get(
                    f"https://rest.uniprot.org/uniprotkb/{entry_id}",
                    params={"format": "json"},
                    timeout=30
                ),
                verbose=verbose
            )

            with open(json_file_path, "w") as outfile:
                json.dump(json_response.json(), outfile, indent=2)
        
        return entry_id
        
    except Exception as e:
        if verbose:
            print(f"Error: {str(e)}")
        return f"Error: {str(e)}"


def score_protein_existence(uniprot_data: Dict[str, Any],
                           verbose: bool = False) -> float:
    """
    Score protein existence evidence (max 15 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for protein existence (0-15)
    """
    protein_existence = uniprot_data.get("proteinExistence", "")
    
    if verbose:
        print(f"\nProtein existence: {protein_existence}")
    
    if protein_existence == "1: Evidence at protein level":
        return 15
    elif protein_existence == "2: Evidence at transcript level":
        return 10
    elif protein_existence == "3: Inferred from homology":
        return 5
    else:  # "4: Predicted" and "5: Uncertain"
        return 0


def score_molecular_function(uniprot_data: Dict[str, Any],
                            verbose: bool = False) -> float:
    """
    Score molecular function annotations (max 20 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for molecular function (0-20)
    """
    molecular_function_terms = set()

    if "uniProtKBCrossReferences" in uniprot_data:
        for ref in uniprot_data["uniProtKBCrossReferences"]:
            if ref.get("database") == "GO" and "properties" in ref:
                for prop in ref["properties"]:
                    if prop.get("key") == "GoTerm" and prop.get("value", "").startswith("F:"):
                        go_term = prop.get("value", "").lower()
                        molecular_function_terms.add(go_term)

    count = len(molecular_function_terms)

    if verbose and molecular_function_terms:
        print(f"GO molecular function terms ({count}):")
        for term in list(molecular_function_terms):
            print(f"- {term}")

    if count >= 4:
        return 20
    elif count == 3:
        return 15
    elif count == 2:
        return 10
    elif count == 1:
        return 5
    else:
        return 0


def score_biological_processes(uniprot_data: Dict[str, Any],
                              verbose: bool = False) -> float:
    """
    Score biological process annotations (max 20 points).

    Evaluates Gene Ontology biological process terms and pathway database
    cross-references (e.g., Reactome) to assess protein's role in biological
    processes.

    Args:
        uniprot_data: UniProt JSON data dictionary
        verbose: If True, prints detailed GO terms and pathways found

    Returns:
        Score for biological processes in range [0, 20]
        - 20 points: 4 or more terms/pathways
        - 15 points: 3 terms/pathways
        - 10 points: 2 terms/pathways
        - 5 points: 1 term/pathway
        - 0 points: No terms/pathways
    """
    biological_process_terms = set()

    if "uniProtKBCrossReferences" in uniprot_data:
        for ref in uniprot_data["uniProtKBCrossReferences"]:
            if ref.get("database") == "GO" and "properties" in ref:
                for prop in ref["properties"]:
                    if prop.get("key") == "GoTerm" and prop.get("value", "").startswith("P:"):
                        go_term = prop.get("value", "").lower()
                        biological_process_terms.add(go_term)
            
            # Also count pathway references
            if ref.get("database") == "Reactome":
                pathway_id = ref.get("id", "")
                biological_process_terms.add(f"pathway:{ref.get('database')}:{pathway_id}")

    count = len(biological_process_terms)

    if verbose and biological_process_terms:
        print(f"GO biological process terms and pathways ({count}):")
        for term in list(biological_process_terms):
            print(f"- {term}")

    if count >= 4:
        return 20
    elif count == 3:
        return 15
    elif count == 2:
        return 10
    elif count == 1:
        return 5
    else:
        return 0


def score_cellular_component(uniprot_data: Dict[str, Any],
                            verbose: bool = False) -> float:
    """
    Score cellular component annotations (max 5 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for cellular component (0-5)
    """
    cellular_component_terms = set()

    if "uniProtKBCrossReferences" in uniprot_data:
        for ref in uniprot_data["uniProtKBCrossReferences"]:
            if ref.get("database") == "GO" and "properties" in ref:
                for prop in ref["properties"]:
                    if prop.get("key") == "GoTerm" and prop.get("value", "").startswith("C:"):
                        go_term = prop.get("value", "").lower()
                        cellular_component_terms.add(go_term)

    count = len(cellular_component_terms)

    if verbose and cellular_component_terms:
        print(f"GO cellular component terms ({count}):")
        for term in list(cellular_component_terms):
            print(f"- {term}")

    if count >= 3:
        return 5
    elif count == 2:
        return 3.5
    elif count == 1:
        return 2
    else:
        return 0


def score_structural_annotation(uniprot_data: Dict[str, Any],
                               verbose: bool = False) -> float:
    """
    Score structural functional annotations (max 10 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for structural annotations (0-10)
    """
    # Define target feature types to count
    target_feature_types = {
        "topological_domain",
        "transmembrane",
        "intramembrane",
        "domain",
        "repeat",
        "zinc_finger",
        "dna_binding",
        "region",
        "coiled_coil",
        "motif",
        "compositional_bias",
        "active_site", 
        "binding_site",
        "site"
    }
    
    # Initialize counters
    feature_type_counts = {feature_type: 0 for feature_type in target_feature_types}
    
    # Count occurrences of each feature type
    if "features" in uniprot_data:
        for feature in uniprot_data["features"]:
            feature_type = feature.get("type", "").lower().replace(" ", "_")
            if feature_type in feature_type_counts:
                feature_type_counts[feature_type] += 1
    
    # Calculate total count of structural annotations
    count = sum(feature_type_counts.values())
    
    # Add detailed counts to verbose output
    if verbose:
        detailed_counts = {k: v for k, v in feature_type_counts.items() if v > 0}
        if detailed_counts:
            print("\nDetailed structural annotation counts:")
            for feature_type, count_detail in detailed_counts.items():
                print(f"- {feature_type.replace('_', ' ').title()}: {count_detail}")
            print(f"Total structural annotations: {count}")
    
    # Assign score based on total count
    if count >= 7:
        return 10
    elif 5 <= count <= 6:
        return 7.5
    elif 3 <= count <= 4:
        return 5
    elif 1 <= count <= 2:
        return 2.5
    else:
        return 0


def score_protein_interactions(uniprot_data: Dict[str, Any],
                              verbose: bool = False) -> float:
    """
    Score protein-protein interactions (max 5 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for protein interactions (0 or 5)
    """
    if "uniProtKBCrossReferences" in uniprot_data:
        for ref in uniprot_data["uniProtKBCrossReferences"]:
            if ref.get("database") in ["IntAct", "BioGRID", "STRING", "ComplexPortal"]:
                if verbose:
                    print(f"Found protein-protein interaction database: {ref.get('database')}")
                return 5
    
    return 0


def score_post_translational_modifications(uniprot_data: Dict[str, Any],
                                          verbose: bool = False) -> float:
    """
    Score post-translational modifications (max 5 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for PTMs (0 or 5)
    """
    ptm_annotations = []

    if "keywords" in uniprot_data:
        for keyword in uniprot_data["keywords"]:
            keyword_category = keyword.get("category", "").strip()
            keyword_name = keyword.get("name", "").strip()
            
            if keyword_category == "PTM" and keyword_name:
                ptm_annotations.append(keyword_name)
                if verbose:
                    print(f"Found PTM keyword: {keyword_name}")

    if ptm_annotations:
        if verbose:
            print(f"Total PTM annotations found: {len(ptm_annotations)}")
        return 5
    else:
        return 0


def score_3d_structures(uniprot_data: Dict[str, Any],
                       verbose: bool = False) -> float:
    """
    Score 3D structure availability (max 5 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for 3D structures (0 or 5)
    """
    if "uniProtKBCrossReferences" in uniprot_data:
        for ref in uniprot_data["uniProtKBCrossReferences"]:
            if ref.get("database") == "PDB":
                if verbose:
                    print("Found 3D structure reference (PDB)")
                return 5
    
    return 0


def score_drug_target(uniprot_data: Dict[str, Any],
                     verbose: bool = False) -> float:
    """
    Score drug target interactions (max 2.5 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for drug targets (0 or 2.5)
    """
    if "uniProtKBCrossReferences" in uniprot_data:
        for ref in uniprot_data["uniProtKBCrossReferences"]:
            if ref.get("database") in ["DrugBank", "ChEMBL", "BindingDB"]:
                if verbose:
                    print(f"Found drug target database: {ref.get('database')}")
                return 2.5
    
    return 0


def score_publication_references(uniprot_data: Dict[str, Any],
                                verbose: bool = False) -> float:
    """
    Score publication references (max 12.5 points).
    
    Parameters
    ----------
    uniprot_data : dict
        UniProt JSON data
    verbose : bool, optional
        Print detailed information
        
    Returns
    -------
    float
        Score for publication references (0-12.5)
    """
    citation_count = len(uniprot_data.get("references", []))
    
    if verbose:
        print(f"Publication references count: {citation_count}")
    
    if citation_count >= 8:
        return 12.5
    elif 6 <= citation_count <= 7:
        return 10
    elif 4 <= citation_count <= 5:
        return 7.5
    elif 2 <= citation_count <= 3:
        return 5
    else:
        return 0


def calculate_protein_score(json_file_path: str,
                           verbose: bool = False) -> float:
    """
    Calculate functional importance score for a protein based on UniProt JSON data.
    
    Parameters
    ----------
    json_file_path : str
        Path to the UniProt JSON file
    verbose : bool, optional
        If True, prints detailed scoring information. Default is False.
        
    Returns
    -------
    float
        Normalized score (0-1) representing protein importance
    """
    # Check if file exists
    if not os.path.exists(json_file_path):
        if verbose:
            print(f"File not found: {json_file_path}")
        return 0.0
    
    try:
        # Load JSON data
        with open(json_file_path, 'r') as f:
            uniprot_data = json.load(f)
        
        # Calculate individual scores using modular functions
        scores = {
            "protein_existence": score_protein_existence(uniprot_data, verbose),
            "molecular_function": score_molecular_function(uniprot_data, verbose),
            "biological_processes": score_biological_processes(uniprot_data, verbose),
            "cellular_component": score_cellular_component(uniprot_data, verbose),
            "structural_annotation": score_structural_annotation(uniprot_data, verbose),
            "protein_interactions": score_protein_interactions(uniprot_data, verbose),
            "post_translational_mods": score_post_translational_modifications(uniprot_data, verbose),
            "3d_structures": score_3d_structures(uniprot_data, verbose),
            "drug_target": score_drug_target(uniprot_data, verbose),
            "publication_references": score_publication_references(uniprot_data, verbose)
        }
        
        # Calculate total score (maximum 100)
        total_score = sum(scores.values())

        # Normalize score (0-1 scale)
        normalized_score = round(total_score / 100, 3)
        
        # Print detailed results if verbose
        if verbose:
            print(f"\nProtein functional importance score: {total_score:.1f}/100")
            print("\nScore breakdown by category:")
            print(f"- Protein existence: {scores['protein_existence']}/15")
            print(f"- Molecular function: {scores['molecular_function']}/20")
            print(f"- Biological processes: {scores['biological_processes']}/20")
            print(f"- Cellular component: {scores['cellular_component']}/5")
            print(f"- Structural annotation: {scores['structural_annotation']}/10")
            print(f"- Protein-protein interactions: {scores['protein_interactions']}/5")
            print(f"- Post-translational modifications: {scores['post_translational_mods']}/5")
            print(f"- 3D structures: {scores['3d_structures']}/5")
            print(f"- Drug target: {scores['drug_target']}/2.5")
            print(f"- Publication references: {scores['publication_references']}/12.5")
        
        # Return only the normalized score (0-1)
        return normalized_score
    
    except Exception as e:
        if verbose:
            print(f"Error calculating protein score: {str(e)}")
        return 0.0


def get_protein_score(taxon_id: str,
                     gene_name: str,
                     verbose: bool = False) -> float:
    """
    Calculate functional importance score for a gene.

    This is the main entry point for protein scoring. It handles the complete
    workflow: fetching UniProt data, caching it, and computing the comprehensive
    functional importance score.

    Args:
        taxon_id: NCBI Taxonomy ID of the organism (e.g., "9606" for human)
        gene_name: Name of the gene to query (e.g., "TP53", "BRCA1")
        verbose: If True, prints detailed progress and scoring breakdown

    Returns:
        Normalized functional importance score in [0, 1] range:
        - 1.0: Highly characterized, functionally important protein
        - 0.5-0.8: Moderately characterized protein
        - 0.0-0.5: Poorly characterized or hypothetical protein
        - 0.0: No UniProt data available

    Notes:
        - Automatically handles gene name conversions (ORF1ab, pol)
        - Results are cached locally to avoid repeated API calls
        - Returns 0.0 if UniProt data cannot be retrieved
    """
    if verbose:
        print(f"\nProcessing taxon ID: {taxon_id}")
    
    try:
        # Retrieve UniProt entry ID
        entry_id = fetch_uniprot_data(taxon_id, gene_name, verbose=verbose)
        
        # Skip if no data found
        if isinstance(entry_id, str) and (entry_id.startswith("No") or entry_id.startswith("Error")):
            if verbose:
                print(f"No UniProt data found for taxon_id={taxon_id}, gene={gene_name}")
            return 0.0

        # Locate the cached JSON file
        current_dir = os.path.dirname(os.path.abspath(__file__))
        json_file_path = os.path.join(current_dir, 'uniprot', f'{entry_id}.json')

        # Calculate the protein score from cached data
        if os.path.exists(json_file_path):
            return calculate_protein_score(json_file_path, verbose=verbose)
        else:
            if verbose:
                print(f"JSON file not found: {json_file_path}")
            return 0.0

    except Exception as e:
        if verbose:
            print(f"Error in get_protein_score: {str(e)}")
        return 0.0


def get_detailed_scores(json_file_path: str,
                       verbose: bool = False) -> Optional[Dict[str, float]]:
    """
    Get detailed breakdown of all scoring components.

    This utility function provides granular access to individual scoring
    components for analysis or debugging purposes.

    Args:
        json_file_path: Path to the cached UniProt JSON file
        verbose: If True, prints detailed information for each component

    Returns:
        Dictionary containing individual scores for each component:
            - protein_existence: 0-15
            - molecular_function: 0-20
            - biological_processes: 0-20
            - cellular_component: 0-5
            - structural_annotation: 0-10
            - protein_interactions: 0-5
            - post_translational_mods: 0-5
            - 3d_structures: 0-5
            - drug_target: 0-2.5
            - publication_references: 0-12.5
            - total_score: 0-100
            - normalized_score: 0-1
        Returns None if file not found or error occurs

    Notes:
        - Useful for understanding which factors contribute to protein score
        - Total score is sum of all components (max 100)
        - Normalized score is total/100 (range [0, 1])
    """
    if not os.path.exists(json_file_path):
        return None
    
    try:
        with open(json_file_path, 'r') as f:
            uniprot_data = json.load(f)
        
        detailed_scores = {
            "protein_existence": score_protein_existence(uniprot_data, verbose),
            "molecular_function": score_molecular_function(uniprot_data, verbose),
            "biological_processes": score_biological_processes(uniprot_data, verbose),
            "cellular_component": score_cellular_component(uniprot_data, verbose),
            "structural_annotation": score_structural_annotation(uniprot_data, verbose),
            "protein_interactions": score_protein_interactions(uniprot_data, verbose),
            "post_translational_mods": score_post_translational_modifications(uniprot_data, verbose),
            "3d_structures": score_3d_structures(uniprot_data, verbose),
            "drug_target": score_drug_target(uniprot_data, verbose),
            "publication_references": score_publication_references(uniprot_data, verbose)
        }
        
        detailed_scores["total_score"] = sum(detailed_scores.values())
        detailed_scores["normalized_score"] = detailed_scores["total_score"] / 100.0
        
        return detailed_scores
        
    except Exception as e:
        if verbose:
            print(f"Error getting detailed scores: {str(e)}")
        return None