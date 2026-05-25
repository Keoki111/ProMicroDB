# DIAMOND-UHGP Annotation Pipeline

[![Python](https://img.shields.io/badge/Python-3.12-blue)](https://www.python.org/) [![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

A lightweight, Python-based pipeline for high-throughput protein sequence alignment using DIAMOND. This tool batch-processes FASTA files against protein databases (e.g., UHGP) and seamlessly merges the alignment results with GTDB taxonomic metadata.

## Key Features

- **Batch Processing**: Automatically scans and processes multiple FASTA/protein sequence files from a designated input directory.

- **Precision Filtering**: Built-in threshold filtering for Sequence Identity, E-value, and Alignment Coverage.

- **Taxonomy Integration**: Automatically joins alignment hits with GTDB metadata. Includes a lineage cleaner to parse standard GTDB strings into discrete taxonomic columns (Domain to Species).

- **Flexible Output Modes**:
  - **Raw**: Pure DIAMOND alignment statistics (fastest, lowest memory).
  - **Standard**: Comprehensive output combining alignment metrics with taxonomy.
  - **Slim**: Minimalist output retaining only Sequence IDs and taxonomic classifications.

## Prerequisites

- Python 3.12
- pandas
- DIAMOND 2.1.17 installed and available in your system `$PATH`.

## Quick Start

1. **Configure**: Open `config.py` to define your database paths (`DIAMOND_DB`, `TAX_FILE`), hardware resources (`threads`), and alignment thresholds.

2. **Input**: Place your query sequence files (e.g., `.fa`, `.fasta`, `.faa`, `.pep`, `.txt`) into the designated `INPUT_FOLDER`.

3. **Prepare**: Ensure the output directory exists (the script does not auto-create it).

4. **Run**: Execute the main script from your terminal:

```bash
python main.py
```

5. **Output**: The filtered and annotated results will be aggregated into a single `.tsv` file in your specified output directory.
