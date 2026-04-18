<div align="center">

<img src="logo.png" alt="RAPID Logo" width="200"/>

# RAPID

### Reproducible Analysis Pipeline for Intensity-based proteomics Data

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Run with Docker](https://img.shields.io/badge/run%20with-docker-0db7ed?logo=docker)](https://www.docker.com/)
[![Run with Conda](https://img.shields.io/badge/run%20with-conda-342B029.svg?logo=anaconda)](https://docs.conda.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI: Lint](https://github.com/ShaunakRaole/RAPID/actions/workflows/lint.yml/badge.svg)](https://github.com/ShaunakRaole/RAPID/actions/workflows/lint.yml)
[![CI: Test](https://github.com/ShaunakRaole/RAPID/actions/workflows/test.yml/badge.svg)](https://github.com/ShaunakRaole/RAPID/actions/workflows/test.yml)

**A fully containerized, reproducible Nextflow pipeline for label-free quantitative proteomics differential abundance analysis.**

[Overview](#overview) •
[Installation](#installation) •
[Quick Start](#quick-start) •
[Pipeline](#pipeline) •
[Parameters](#parameters) •
[Output](#output) •
[Test Data](#test-data) •
[Citation](#citation)

</div>

---

## Overview

RAPID is a reproducible bioinformatics pipeline for analyzing label-free quantitative (LFQ) proteomics data from MaxQuant output. It implements the statistically rigorous workflow recommended by [Peng et al. (2024)](https://doi.org/10.1038/s41467-024-47899-w), which demonstrated that tool choice and normalization strategy substantially affect differential abundance results in proteomics.

RAPID addresses the reproducibility crisis in proteomics data analysis by providing:

- **Single-command execution** from raw MaxQuant output to publication-ready results
- **Full containerization** via Docker — identical results on any system
- **Statistically appropriate methods** — DEqMS accounts for peptide-level quantification uncertainty, unlike generic limma
- **Multi-organism support** — human, mouse, rat, and yeast pathway enrichment via clusterProfiler
- **Automated validation** — built-in output integrity checks on every run
- **CI/CD integration** — GitHub Actions lint and end-to-end tests on every commit

### Why RAPID?

| Tool | Reproducible | Containerized | Scriptable | DEqMS | Multi-organism enrichment |
|------|-------------|---------------|------------|-------|--------------------------|
| Perseus | ❌ | ❌ | ❌ | ❌ | ❌ |
| Custom R script | ⚠️ | ❌ | ✅ | ⚠️ | ⚠️ |
| **RAPID** | ✅ | ✅ | ✅ | ✅ | ✅ |

---

## Pipeline

RAPID implements an 8-module workflow:

```
proteinGroups.txt
      │
      ▼
┌─────────────────┐
│  1. DOWNLOAD    │  wget from URL (PRIDE, GitHub, etc.)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  2. QC_FILTER   │  Remove contaminants, reverse hits, low-coverage proteins
└────────┬────────┘
         │
    ┌────┴─────┐
    ▼          ▼
┌───────┐  ┌──────────┐
│ 3.    │  │  4. PCA  │  Quality control visualization
│NORMAL │  └────┬─────┘
│  IZE  │       │
└───┬───┘       │
    │            │
    ▼            │
┌──────────┐     │
│ 5.DA_    │     │
│ MULTI    │  All pairwise DEqMS comparisons
└────┬─────┘     │
     │            │
     ▼            │
┌──────────┐      │
│ 6. ENRI- │      │
│  CHMENT  │  GO/KEGG pathway enrichment
└────┬─────┘      │
     │             │
     └──────┬──────┘
            ▼
     ┌─────────────┐
     │  7. REPORT  │  Self-contained HTML report
     └──────┬──────┘
            │
            ▼
     ┌─────────────┐
     │  8. VALIDATE│  Output integrity checks
     └─────────────┘
```

### Modules

| # | Module | Script | Description |
|---|--------|--------|-------------|
| 1 | `DOWNLOAD_DATA` | — | Downloads `proteinGroups.txt` from any URL |
| 2 | `QC_FILTER` | `01_qc_filter.R` | Removes contaminants (`+`), reverse hits (`+`), site-only proteins (`+`), and proteins below minimum valid value threshold |
| 3 | `NORMALIZE` | `02_normalize.R` | log2 transformation followed by median centering; generates before/after boxplot |
| 4 | `PCA` | `03_pca.R` | PCA of normalized matrix with row-mean imputation for missing values |
| 5 | `DA_MULTI` | `04_da_multi.R` | All pairwise DEqMS differential abundance comparisons; one volcano plot per comparison |
| 6 | `ENRICHMENT` | `05_enrichment.R` | GO Biological Process and KEGG enrichment via clusterProfiler; organism-aware |
| 7 | `REPORT` | `06_report.R` | Self-contained HTML report with all figures, tables, and reproducibility info |
| 8 | `VALIDATE` | `08_validate.R` | Automated output validation — checks file existence, dimensions, column integrity, p-value ranges, and report structure |

---

## Installation

### Requirements

| Dependency | Version | Purpose |
|------------|---------|---------|
| [Nextflow](https://www.nextflow.io/) | ≥ 23.04.0 | Workflow management |
| [Docker](https://www.docker.com/) | any | Recommended execution environment |
| [Conda](https://docs.conda.io/) | any | Alternative execution environment |
| Java | ≥ 11 | Required by Nextflow |

### Step 1 — Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

### Step 2 — Clone the repository

```bash
git clone https://github.com/ShaunakRaole/RAPID.git
cd RAPID
```

### Step 3 — Set up your environment

**Option A: Docker (recommended)**

Pull the pre-built image:

```bash
docker pull shaunakraole/rapid:1.0
```

Or build locally:

```bash
docker build -t shaunakraole/rapid:1.0 .
```

**Option B: Conda**

```bash
conda env create -f environment.yml
conda activate rapid
```

> **Note:** `org.Mm.eg.db`, `org.Rn.eg.db`, and `org.Sc.sgd.db` (mouse, rat, yeast annotation packages) must be installed manually after conda env creation due to a known Bioconductor post-install issue on macOS ARM64:
> ```bash
> conda activate rapid
> Rscript -e "BiocManager::install(c('org.Mm.eg.db', 'org.Rn.eg.db', 'org.Sc.sgd.db'))"
> ```
> These packages are pre-installed in the Docker image and are not required if using Docker.

---

## Quick Start

### Run with Docker (recommended)

```bash
nextflow run main.nf -profile docker \
    --data_url    "https://raw.githubusercontent.com/statOmics/PDA21/data/quantification/cptacAvsB_lab3/proteinGroups.txt" \
    --sample_sheet test-data/test_sample_sheet.csv \
    --condition_a  6A \
    --condition_b  6B \
    --organism     human \
    --outdir       results
```

### Run with Conda

```bash
nextflow run main.nf -profile conda \
    --data_url    "https://your-pride-url/proteinGroups.txt" \
    --sample_sheet my_sample_sheet.csv \
    --condition_a  control \
    --condition_b  treatment \
    --organism     human \
    --outdir       results
```

### Run the bundled test profile

```bash
nextflow run main.nf -profile test,docker
```

---

## Sample Sheet

The sample sheet is a plain CSV file describing your experimental design. It must have exactly three columns:

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | ✅ | Unique sample identifier used in output files |
| `condition` | ✅ | Experimental condition — must match `--condition_a` / `--condition_b` |
| `file_name` | ✅ | Must exactly match the suffix after `LFQ intensity ` in the MaxQuant `proteinGroups.txt` column header |

**Example:**

```csv
sample,condition,file_name
control_rep1,control,ctrl_1
control_rep2,control,ctrl_2
control_rep3,control,ctrl_3
treatment_rep1,treatment,treat_1
treatment_rep2,treatment,treat_2
treatment_rep3,treatment,treat_3
```

**How to find your `file_name` values:**

```bash
# Print all LFQ column names from your proteinGroups.txt
head -1 proteinGroups.txt | tr '\t' '\n' | grep "LFQ intensity"
# The file_name is everything after "LFQ intensity "
```

> **Multi-condition support:** RAPID automatically detects all conditions in the sample sheet and runs all pairwise comparisons. You do not need to specify `--condition_a` and `--condition_b` for 3+ condition datasets — just leave them unset.

---

## Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--data_url` | `string` | URL to MaxQuant `proteinGroups.txt` file. Must be publicly accessible (PRIDE FTP, GitHub raw, etc.) |
| `--sample_sheet` | `path` | Path to sample sheet CSV file (see [Sample Sheet](#sample-sheet) format) |

### Analysis Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--condition_a` | `string` | `null` | Reference condition name for labelling. Optional for 3+ condition datasets |
| `--condition_b` | `string` | `null` | Comparison condition name for labelling. Optional for 3+ condition datasets |
| `--min_valid_values` | `integer` | `2` | Minimum number of samples with a valid (non-zero, non-NA) LFQ intensity for a protein to be retained after QC filtering |
| `--fdr_cutoff` | `float` | `0.05` | Benjamini-Hochberg adjusted p-value threshold for significance. Applied to both DA results and pathway enrichment |
| `--fc_cutoff` | `float` | `1.0` | log2 fold change threshold for volcano plot and significance classification. Proteins must exceed both `--fdr_cutoff` AND `--fc_cutoff` to be considered significant |
| `--organism` | `string` | `human` | Organism for pathway enrichment annotation. Options: `human`, `mouse`, `rat`, `yeast` |

### Infrastructure Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--outdir` | `string` | `results` | Directory for all pipeline outputs. Will be created if it does not exist |

### Execution Profiles

| Profile | Description | Use case |
|---------|-------------|----------|
| `docker` | Uses `shaunakraole/rapid:1.0` container | **Recommended** — fully reproducible, works on any system with Docker |
| `conda` | Uses `environment.yml` conda environment | Local development on Linux/macOS |
| `test` | Pre-configured test run using CPTAC benchmark data | Verifying installation, CI/CD |
| `test,docker` | Test profile + Docker | Clean end-to-end validation |

---

## Output

All outputs are written to `--outdir` (default: `results/`):

```
results/
├── raw/
│   └── proteinGroups.txt              # Downloaded raw MaxQuant output
│
├── qc/
│   ├── filtered_matrix.tsv            # QC-filtered intensity matrix (tab-separated)
│   ├── qc_summary.txt                 # QC statistics (proteins removed, retained)
│   ├── normalization_boxplot.png      # Intensity distributions before/after normalization
│   ├── pca_plot.png                   # PCA of normalized samples coloured by condition
│   └── pca_coords.tsv                 # PCA coordinates per sample
│
├── normalized/
│   └── normalized_matrix.tsv          # log2 + median-centered intensity matrix
│
├── results/
│   ├── da_results_all_comparisons.tsv # Combined DA results across all pairwise comparisons
│   ├── da_results_<B>_vs_<A>.tsv      # Per-comparison DA results
│   └── volcano_<B>_vs_<A>.png         # Volcano plot per comparison
│
├── enrichment/
│   ├── enrichment_results.tsv         # Combined GO and KEGG enrichment results
│   ├── go_dotplot_<comparison>.png    # GO Biological Process dot plot (if significant)
│   └── kegg_dotplot_<comparison>.png  # KEGG pathway dot plot (if significant)
│
├── report/
│   └── proteomics_report.html         # Self-contained HTML report (all figures embedded)
│
├── validation/
│   └── validation_report.txt          # Automated output validation results
│
├── pipeline_report.html               # Nextflow execution report
├── timeline.html                      # Process execution timeline
└── pipeline_dag.html                  # Pipeline DAG visualization
```

### DA Results Columns

The main output `da_results_all_comparisons.tsv` contains:

| Column | Description |
|--------|-------------|
| `protein_id` | UniProt protein identifier(s) |
| `gene_names` | Gene symbol(s) |
| `log2FC` | log2 fold change (condition_b / condition_a) |
| `pvalue` | DEqMS sca.P.Value |
| `adj_pvalue` | BH-adjusted p-value |
| `significant` | `TRUE` if passes both FDR and FC thresholds |
| `direction` | `Up`, `Down`, or `NS` |
| `comparison` | Comparison label (e.g. `treatment_vs_control`) |
| `condition_a` | Reference condition |
| `condition_b` | Comparison condition |

---

## Test Data

RAPID ships with a pre-configured test profile using the CPTAC benchmark dataset — a well-characterised proteomics gold standard used extensively in methods development.

**Dataset:** CPTAC Study 6 — UPS1 spike-in (48 human proteins) in *S. cerevisiae* background  
**Source:** [statOmics/PDA21](https://github.com/statOmics/PDA21) — Goeminne et al. (2016)  
**Conditions:** 6A (0.25 fmol/μL spike-in) vs 6B (0.74 fmol/μL spike-in)  
**Samples:** 6 total (3 replicates per condition)  
**Expected runtime:** ~5 minutes on a laptop  
**Expected output:** ~2 significant UPS1 human proteins; enrichment sparse (mixed yeast/human background)

```bash
# Run the test profile
nextflow run main.nf -profile test,docker

# Expected outputs in test_results/
```

**Test sample sheet** is available at `test-data/test_sample_sheet.csv`.

---

## Scientific Background

### Why DEqMS instead of limma?

Standard limma assumes equal prior variance across all proteins. In proteomics, proteins quantified by more peptides or PSMs have inherently more reliable intensity estimates. DEqMS (Zhen et al., 2020) explicitly models this relationship, providing better statistical power and accuracy — particularly important for proteins quantified by a single peptide.

### Why median centering normalization?

Median centering is a robust, assumption-free normalization method appropriate for LFQ proteomics data where systematic shifts between samples are expected but large-scale biological differences may exist. It preserves relative differences between samples better than quantile normalization, which can introduce artifacts when conditions have genuinely different proteome compositions.

### Reproducibility

Tool choice in proteomics differential abundance analysis substantially affects results. Peng et al. (2024) showed that different pipelines applied to the same data can produce drastically different lists of significant proteins. RAPID implements the recommended workflow with fully pinned software versions, ensuring that results are reproducible across systems, time, and users.

---

## Dependencies

### R Packages

| Package | Version | Source | Purpose |
|---------|---------|--------|---------|
| `data.table` | ≥1.14 | CRAN | Fast data I/O and manipulation |
| `ggplot2` | ≥3.4 | CRAN | Visualization |
| `ggrepel` | ≥0.9 | CRAN | Non-overlapping labels on volcano plots |
| `rmarkdown` | ≥2.20 | CRAN | HTML report rendering |
| `knitr` | ≥1.40 | CRAN | Report code execution |
| `kableExtra` | ≥1.3 | CRAN | HTML table formatting |
| `bit64` | ≥4.0 | CRAN | 64-bit integer support for large intensity values |
| `limma` | ≥3.54 | Bioconductor | Linear model framework for DA analysis |
| `DEqMS` | ≥1.16 | Bioconductor | Peptide-count-aware variance estimation |
| `clusterProfiler` | ≥4.6 | Bioconductor | GO and KEGG enrichment analysis |
| `enrichplot` | ≥1.18 | Bioconductor | Enrichment result visualization |
| `org.Hs.eg.db` | ≥3.16 | Bioconductor | Human gene annotation |
| `org.Mm.eg.db` | ≥3.16 | Bioconductor | Mouse gene annotation |
| `org.Rn.eg.db` | ≥3.16 | Bioconductor | Rat gene annotation |
| `org.Sc.sgd.db` | ≥3.16 | Bioconductor | Yeast gene annotation |

---

## Troubleshooting

### Common issues

**`LFQ column mismatch` error in QC_FILTER**
> The `file_name` column in your sample sheet does not match the column names in `proteinGroups.txt`. Run `head -1 proteinGroups.txt | tr '\t' '\n' | grep "LFQ"` to see exact column names and update your sample sheet accordingly.

**`make.names()` condition name warning**
> Condition names starting with numbers (e.g. `6A`, `6B`) are automatically sanitized. This is expected behaviour — results will be correct.

**Empty enrichment results**
> This is expected when: (a) the dataset contains non-human proteins (e.g. yeast), (b) fewer than 3 significant proteins are detected, or (c) the significant proteins do not cluster into annotated pathways at the chosen FDR threshold. Try `--fdr_cutoff 0.1 --fc_cutoff 0.5` for more permissive thresholds.

**`org.Mm.eg.db` fails to install via conda**
> Known Bioconductor post-link script issue on macOS ARM64. Install manually: `Rscript -e "BiocManager::install('org.Mm.eg.db')"`. All packages are pre-installed in the Docker image.

**Pipeline runs but report is empty**
> Delete the `work/` directory and rerun — stale cached results from a previous failed run may be staged. `rm -rf work/ && nextflow run main.nf ...`

---

## Repository Structure

```
RAPID/
├── .github/
│   └── workflows/
│       ├── lint.yml          # R syntax, Nextflow DSL, Dockerfile linting
│       └── test.yml          # End-to-end pipeline test with output validation
│
├── bin/
│   ├── 01_qc_filter.R        # QC and contaminant filtering
│   ├── 02_normalize.R        # log2 transformation and median centering
│   ├── 03_pca.R              # PCA of normalized matrix
│   ├── 04_da_multi.R         # All pairwise DEqMS differential abundance
│   ├── 05_enrichment.R       # GO/KEGG pathway enrichment (organism-aware)
│   ├── 06_report.R           # HTML report generation
│   └── 08_validate.R         # Output validation
│
├── test-data/
│   └── test_sample_sheet.csv # Sample sheet for CPTAC test dataset
│
├── .gitignore
├── Dockerfile                # Docker image definition (rocker/r-ver:4.3.3 base)
├── LICENSE
├── README.md
├── environment.yml           # Conda environment specification
├── logo.png                  # RAPID pipeline logo
├── main.nf                   # Nextflow DSL2 workflow
└── nextflow.config           # Pipeline configuration and profiles
```

---

## Contributing

Contributions are welcome. Please open an issue before submitting a pull request to discuss proposed changes.

1. Fork the repository
2. Create a feature branch (`git checkout -b feat/my-feature`)
3. Make your changes with atomic commits and descriptive messages
4. Ensure all CI checks pass
5. Open a pull request against `main`

---

## Citation

If you use RAPID in your research, please cite:

**RAPID Pipeline**
> Raole, S. (2026). RAPID: Reproducible Analysis Pipeline for Intensity-based proteomics Data. Georgia Institute of Technology. https://github.com/ShaunakRaole/RAPID

**DEqMS**
> Zhen, Y., Aardema, M.L., Medema, E.M., Lodder, E.M., Linke, W.A., de Weger, R.A., & Chamuleau, M.E.D. (2020). DEqMS: A Method for Accurate Variance Estimation in Differential Protein Expression Analysis. *Molecular & Cellular Proteomics*, 19(6), 1047–1057. https://doi.org/10.1074/mcp.TIR119.001646

**clusterProfiler**
> Yu, G., Wang, L.G., Han, Y., & He, Q.Y. (2012). clusterProfiler: an R Package for Comparing Biological Themes Among Gene Clusters. *OMICS: A Journal of Integrative Biology*, 16(5), 284–287. https://doi.org/10.4253/omics.20120018

**Nextflow**
> Di Tommaso, P., Chatzou, M., Floden, E.W., Barja, P.P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316–319. https://doi.org/10.1038/nbt.3820

**OpDEA benchmarking**
> Peng, M., Guo, T., & Hermjakob, H. (2024). Optimizing differential expression analysis for proteomics data via high-performing rules and ensemble inference. *Nature Communications*, 15, 3515. https://doi.org/10.1038/s41467-024-47899-w

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

<div align="center">

Developed by **Shaunak Raole** · Georgia Institute of Technology · BIOL 8802

</div>
