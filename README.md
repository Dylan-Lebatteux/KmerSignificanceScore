# K-mer Significance Score (KSS)

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A bioinformatics tool for analyzing viral genomic sequences using k-mer based significance scoring. This method integrates discriminative, mutational, and protein-level features to identify biologically relevant sequence variations.

## 🎯 Overview

The K-mer Significance Score (KSS) framework analyzes viral sequence datasets to identify functionally important genomic positions by combining:

- **Discriminative Score**: Quantifies how well a k-mer distinguishes between viral strains/classes
- **Mutational Score**: Evaluates the biochemical impact of amino acid substitutions
- **Protein Score**: Assesses functional relevance using structural and Gene Ontology annotations

### Key Features

✅ **Incremental Storage** - Stores only mutations (70-99% file size reduction)
✅ **Scalable** - Handles 300k+ sequences efficiently
✅ **Publication-Ready** - VCF-like JSON format with metadata
✅ **Configurable** - YAML configuration per dataset
✅ **Optimized** - O(1) lookups, fast compilation

## 📦 Installation

### Requirements

- Python 3.8 or higher
- Dependencies listed in `requirements.txt`

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/KSS.git
cd KSS

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## 🚀 Quick Start

### Basic Usage

```bash
# Analyze all configured datasets
python main.py

# Analyze specific dataset
python main.py data/Human_betaherpesvirus_5/config.yaml

# Analyze multiple datasets
python main.py data/Human_betaherpesvirus_5/config.yaml data/SARS-CoV-2/config.yaml

# Use glob patterns
python main.py data/*/config.yaml
```

### Configuration

Each dataset has its own `config.yaml` file:

```yaml
# data/Human_betaherpesvirus_5/config.yaml
dataset:
  name: "Human_betaherpesvirus_5"
  genes: ["UL55", "UL73", "US28", "UL33"]

parameters:
  k: 9
  threshold: 25

  weights:
    discriminative: 1
    mutational: 1
    protein: 1

  scoring:
    discriminative_threshold: 0.9
    n_theoretical_discriminable_classes: 4
    mutational_matrix: "MIYATA_EVO"
    substitution_matrix: "BLOSUM62"

  alignment:
    open_gap_score: -10
    extend_gap_score: -1

  output:
    evaluation: true
    save_raw: true
```

### Data Format

Input data should be organized as:

```
data/
└── Virus_Name/
    ├── config.yaml           # Dataset configuration
    ├── GeneA/
    │   ├── reference.gb      # GenBank reference
    │   └── sequences.fasta   # Aligned sequences (class labeled)
    └── GeneB/
        ├── reference.gb
        └── sequences.fasta
```

**Sequence headers must include class labels:** `>sequence_id|Virus|Gene|CLASS`

## 📁 Project Structure

```
KSS/
├── main.py                      # Main analysis script
├── requirements.txt             # Python dependencies
├── LICENSE                      # MIT License
│
├── src/                         # Source code
│   ├── kanalyzer.py            # K-mer analysis engine
│   ├── kss.py                  # KSS score computation
│   ├── discriminative_score.py # Discriminative scoring
│   ├── mutation_score.py       # Mutational impact scoring
│   ├── protein_score.py        # Protein-level scoring
│   ├── utils.py                # Utility functions
│   ├── substitution_matrices/  # AA substitution matrices
│   └── uniprot/                # UniProt data cache
│
├── data/                        # Viral genome datasets
│   ├── Human_betaherpesvirus_5/
│   │   ├── config.yaml
│   │   ├── UL55/, UL73/, US28/, UL33/
│   │   └── results/
│   ├── Human_immunodeficiency_virus_1/
│   │   ├── config.yaml
│   │   ├── gag/, pol/, env/
│   │   └── results/
│   └── Severe_acute_respiratory_syndrome_coronavirus_2/
│       ├── config.yaml
│       ├── ORF1ab/, S/, M/, N/, E/
│       └── results/
│
├── notebooks/                   # Jupyter analysis notebooks
├── examples/                    # Example outputs
└── archive/                     # Archived docs & tests
```

## 📊 Output Files

Results are saved in `data/Virus_Name/results/`:

| File | Description |
|------|-------------|
| `*_results.json` | Raw results with incremental storage (mutations only) |
| `*_compiled_results.json` | Compiled results with KSS scores |
| `*_mi_scores.json` | Mutual information scores (if evaluation=true) |
| `*_chi2_scores.json` | Chi-squared test results (if evaluation=true) |
| `*_odds_ratio_scores.json` | Odds ratio calculations (if evaluation=true) |

### Result Structure

**Raw results** (`*_results.json`):
```json
{
  "_metadata": {
    "version": "2.0.0",
    "tool": "K-mer Significance Score (KSS)",
    "k": 9,
    "sequences_analyzed": 3124
  },
  "genes": {
    "E": {
      "metadata": {
        "reference_length": 228,
        "num_sequences": 3124
      },
      "sequences": {
        "seq_id|Virus|Gene|B.1.1": {
          "class": "B.1.1",
          "mutations": [
            {
              "position": 19,
              "ref_kmer": "CATGCATGC",
              "alt_kmer": "ATGCTTGCA",
              "aa_changes": [
                {"notation": "C7T", "mut_score": 0.85}
              ]
            }
          ]
        }
      }
    }
  }
}
```

**Compiled results** (`*_compiled_results.json`):
```json
{
  "E": {
    "19": {
      "ref": "CATGCATGC",
      "alts": {
        "ATGCTTGCA": {
          "class_counts": {"B.1.1": 45, "B.1.2": 12},
          "amino_acid_changes": {"C7T": 0.85}
        }
      },
      "discriminative_score": 0.82,
      "mutational_score": 0.75,
      "protein_score": 0.65,
      "kss": 0.74
    }
  }
}
```

## 🧬 Datasets

The repository includes GenBank reference genomes and configuration files for three viral systems. **Sequence FASTA files are not included due to their size** (up to 5.6 GB per gene). To obtain the datasets:

1. **SARS-CoV-2** (278,738 sequences, 19 variants) — Download from [GISAID](https://gisaid.org/) or [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/)
2. **HIV-1** (12,223 sequences, 15 subtypes) — Download from [Los Alamos HIV Sequence Database](https://www.hiv.lanl.gov/)
3. **HCMV** (399–646 sequences, 4–8 genotypes) — Download from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)

Place FASTA files as `data/<Virus>/<Gene>/sequences.fasta` with headers formatted as `>sequence_id|Virus|Gene|CLASS`.

## 📈 Performance

### File Size Reduction (Incremental Storage)

| Virus | Sequences | Old Format | New Format | Reduction |
|-------|-----------|------------|------------|-----------|
| SARS-CoV-2 | 3,124 | 1.3 GB | 70 MB | **98.5%** ✅ |
| CMV | 367 | - | - | ~50% |
| HIV-1 | 367 | 2.0 GB | 2.0 GB | Minimal* |

*HIV-1 has very high variability (~90% positions mutated), so minimal reduction is expected.

### Scalability

- **300k sequences:** ~7 GB (manageable!)
- **Typical runtime:** 10-15 minutes for 3 datasets (367 sequences each)
- **Memory usage:** 3-5 GB peak
- **Processing:** Sequential (optimal for I/O-bound workload)

## ⚙️ Advanced Configuration

### Substitution Matrices

Available matrices for mutational scoring:
- `MIYATA_EVO` (default) - Evolutionary distance
- `BLOSUM62` - BLOSUM substitution matrix
- `GRANTHAM` - Chemical distance
- `PAM250` - Point accepted mutations

### Protein Scoring Components

The protein score (0-100 points) integrates:
- Protein existence evidence (15 pts)
- Molecular function GO terms (20 pts)
- Biological process GO terms (20 pts)
- Cellular component GO terms (5 pts)
- Structural annotations (10 pts)
- Protein-protein interactions (5 pts)
- Post-translational modifications (5 pts)
- 3D structure availability (5 pts)
- Drug target information (2.5 pts)
- Publication references (12.5 pts)

## 🔍 Workflow

```
1. Sequence Analysis (kanalyzer.py)
   ├─ Read GenBank reference
   ├─ Read FASTA query sequences
   ├─ Align proteins → back-translate to nucleotides
   └─ Extract k-mers and identify mutations (incremental storage)

2. Result Compilation (kss.py)
   ├─ Aggregate mutations by position/k-mer
   ├─ Count occurrences by class
   ├─ Count reference k-mers (sequences without mutation)
   └─ Filter by threshold (25% in at least one class)

3. KSS Score Calculation (kss.py)
   ├─ Discriminative score (class separation)
   ├─ Mutational score (substitution matrices)
   ├─ Protein score (UniProt annotations)
   └─ Weighted final KSS score

4. Statistical Evaluation (optional)
   ├─ Mutual Information
   ├─ Chi-squared test
   └─ Odds Ratios
```

## 🐛 Troubleshooting

### Common Issues

**Issue:** `UnicodeEncodeError` on Windows
```bash
# Already fixed in main.py with UTF-8 encoding wrapper
```

**Issue:** Out of memory with large datasets
```bash
# Solution: Process one dataset at a time
python main.py data/virus1/config.yaml
python main.py data/virus2/config.yaml
```

**Issue:** Slow performance
```bash
# Check: Sequential processing is optimal (don't use multiprocessing)
# The workload is I/O-bound, not CPU-bound
```

## 📚 Documentation

For detailed documentation:
- See `notebooks/` for analysis examples and methodology validation
- Check `archive/documentation/` for archived technical documentation
- Read `CODE_VERIFICATION.md` for code quality report

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

For questions or issues, please open an issue on GitHub.

## 🙏 Acknowledgments

- BioPython team for sequence analysis tools
- UniProt for protein annotations
- NCBI for genomic data

---

**Version:** 2.0.0
**Last Updated:** 2025-11-08
**Status:** Production-ready 🚀
