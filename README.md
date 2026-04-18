# ProteomicsDA

A reproducible Nextflow pipeline for proteomics differential abundance analysis.

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-342B029.svg?logo=anaconda)](https://docs.conda.io/en/latest/)

## Overview

ProteomicsDA implements a reproducible 5-step workflow for proteomics
differential abundance analysis starting from a MaxQuant `proteinGroups.txt`
file:

| Step | Module | Description |
|------|--------|-------------|
| 1 | `DOWNLOAD_DATA` | Fetch `proteinGroups.txt` from a public URL (e.g. PRIDE) |
| 2 | `QC_FILTER`     | Remove contaminants, reverse hits, low-coverage proteins |
| 3 | `NORMALIZE`     | log2 transform + median centering |
| 4 | `DA_ANALYSIS`   | DEqMS differential abundance analysis |
| 5 | `REPORT`        | Self-contained HTML report with volcano plot and results table |

The pipeline is grounded in recommendations from Peng et al. (2024, Nature
Communications), which showed that tool choice substantially affects proteomics
differential abundance results. ProteomicsDA implements the recommended
workflow for DDA data: directLFQ-style normalization + DEqMS.

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 23.04.0
- [Docker](https://www.docker.com/) (recommended) **or** [Conda](https://docs.conda.io/)

## Quick start

```bash
# 1. Install Nextflow
curl -s https://get.nextflow.io | bash

# 2. Run with Docker (recommended)
nextflow run main.nf \
    -profile docker \
    --data_url    "https://your-pride-url/proteinGroups.txt" \
    --sample_sheet data/your_sample_sheet.csv \
    --condition_a  control \
    --condition_b  treatment

# 3. Run with Conda
nextflow run main.nf \
    -profile conda \
    --data_url    "https://your-pride-url/proteinGroups.txt" \
    --sample_sheet data/your_sample_sheet.csv \
    --condition_a  control \
    --condition_b  treatment

# 4. Run with bundled test data
nextflow run main.nf -profile test,docker
```

## Sample sheet format

The sample sheet is a CSV file with three required columns:

| Column | Description |
|--------|-------------|
| `sample` | Unique sample name used in output files |
| `condition` | Experimental condition (must match `--condition_a` / `--condition_b`) |
| `file_name` | Must match the suffix of the `LFQ intensity <file_name>` column in MaxQuant output |

Example:

```csv
sample,condition,file_name
control_1,control,control_1
control_2,control,control_2
treatment_1,treatment,treatment_1
treatment_2,treatment,treatment_2
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--data_url` | required | URL to MaxQuant proteinGroups.txt |
| `--sample_sheet` | required | Path to sample sheet CSV |
| `--condition_a` | required | Reference condition name |
| `--condition_b` | required | Comparison condition name |
| `--min_valid_values` | `2` | Minimum samples with valid intensity per protein |
| `--fdr_cutoff` | `0.05` | BH-adjusted p-value threshold |
| `--fc_cutoff` | `1.0` | log2 fold change threshold for volcano plot |
| `--outdir` | `results` | Output directory |

## Output

```
results/
├── raw/
│   └── proteinGroups.txt          # Downloaded raw data
├── qc/
│   ├── filtered_matrix.tsv        # QC-filtered intensity matrix
│   └── qc_summary.txt             # QC statistics
├── normalized/
│   ├── normalized_matrix.tsv      # Log2-normalized matrix
│   └── normalization_boxplot.pdf  # Before/after distribution plots
├── results/
│   ├── da_results.tsv             # Full differential abundance results
│   └── volcano_plot.pdf           # Volcano plot
├── report/
│   └── proteomics_report.html     # Self-contained HTML report
├── pipeline_report.html           # Nextflow execution report
├── timeline.html                  # Nextflow timeline
└── pipeline_dag.html              # Pipeline DAG visualization
```

## Building the Docker image

```bash
docker build -t shaunak/proteomicsda:1.0 .
```

## References

- Zhi et al. (2020). DEqMS: A Method for Accurate Variance Estimation in
  Differential Protein Expression Analysis. *Molecular & Cellular Proteomics*.
  doi:10.1074/mcp.TIR119.001646

- Peng et al. (2024). Optimizing differential expression analysis for proteomics
  data via high-performing rules and ensemble inference. *Nature Communications*.
  doi:10.1038/s41467-024-47899-w

- Di Tommaso et al. (2017). Nextflow enables reproducible computational workflows.
  *Nature Biotechnology*. doi:10.1038/nbt.3820
