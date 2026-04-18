#!/usr/bin/env Rscript
# =============================================================================
# 04_report.R
# Generates a self-contained HTML report summarising the pipeline run
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rmarkdown)
  library(knitr)
})

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    qc_summary   = NULL,
    boxplot      = NULL,
    da_results   = NULL,
    volcano      = NULL,
    sample_sheet = NULL,
    condition_a  = NULL,
    condition_b  = NULL,
    fdr_cutoff   = 0.05,
    fc_cutoff    = 1.0,
    output       = "proteomics_report.html"
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--qc-summary"   = { defaults$qc_summary   <- args[i+1]; i <- i+2 },
      "--boxplot"      = { defaults$boxplot       <- args[i+1]; i <- i+2 },
      "--da-results"   = { defaults$da_results    <- args[i+1]; i <- i+2 },
      "--volcano"      = { defaults$volcano       <- args[i+1]; i <- i+2 },
      "--sample-sheet" = { defaults$sample_sheet  <- args[i+1]; i <- i+2 },
      "--condition-a"  = { defaults$condition_a   <- args[i+1]; i <- i+2 },
      "--condition-b"  = { defaults$condition_b   <- args[i+1]; i <- i+2 },
      "--fdr-cutoff"   = { defaults$fdr_cutoff    <- as.numeric(args[i+1]); i <- i+2 },
      "--fc-cutoff"    = { defaults$fc_cutoff     <- as.numeric(args[i+1]); i <- i+2 },
      "--output"       = { defaults$output        <- args[i+1]; i <- i+2 },
      {
        stop(sprintf("Unknown argument: %s", args[i]))
      }
    )
  }

  for (req in c("qc_summary", "da_results", "sample_sheet",
                "condition_a", "condition_b")) {
    if (is.null(defaults[[req]])) {
      stop(sprintf("--%s is required", gsub("_", "-", req)))
    }
  }
  defaults
}

opt <- parse_args(args)

# -----------------------------------------------------------------------------
# Build inline Rmd template and render
# We write a temp .Rmd file and render it with the paths as parameters
# -----------------------------------------------------------------------------

rmd_template <- '---
title: "ProteomicsDA — Analysis Report"
date: "`r format(Sys.time(), \'%B %d, %Y\')`"
output:
  html_document:
    theme: flatly
    toc: true
    toc_float: true
    self_contained: true
    code_folding: hide
params:
  qc_summary:   ""
  boxplot_path: ""
  da_results:   ""
  volcano_path: ""
  sample_sheet: ""
  condition_a:  ""
  condition_b:  ""
  fdr_cutoff:   0.05
  fc_cutoff:    1.0
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(data.table)
library(knitr)
library(kableExtra)
library(ggplot2)
```

## Overview

This report summarises a reproducible proteomics differential abundance analysis
comparing **`r params$condition_b`** versus **`r params$condition_a`**.

The pipeline was run using [ProteomicsDA](https://github.com/sraole3/proteomicsda),
a Nextflow-based workflow implementing:

- **QC filtering**: removal of contaminants, reverse hits, and low-coverage proteins
- **Normalization**: log2 transformation and median centering
- **Differential abundance**: DEqMS (Bioconductor)

---

## 1. Sample Overview

```{r sample-table}
ss <- fread(params$sample_sheet)
kable(ss, caption = "Samples included in this analysis") |>
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)
```

---

## 2. QC Summary

```{r qc-summary}
qc_lines <- readLines(params$qc_summary)
cat(paste(qc_lines, collapse = "\n"))
```

---

## 3. Normalization

The plot below shows protein intensity distributions before and after
log2 transformation and median centering. After normalization, sample
medians should be aligned.

```{r boxplot, fig.cap="Protein intensity distributions before and after normalization", out.width="100%"}
knitr::include_graphics(params$boxplot_path)
```

---

## 4. Differential Abundance Results

Proteins were tested for differential abundance between
**`r params$condition_b`** and **`r params$condition_a`** using DEqMS.
Significance threshold: FDR < `r params$fdr_cutoff`,
|log2FC| ≥ `r params$fc_cutoff`.

### 4.1 Volcano plot

```{r volcano, fig.cap="Volcano plot. Red = upregulated in condition B; Blue = downregulated.", out.width="100%"}
knitr::include_graphics(params$volcano_path)
```

### 4.2 Summary statistics

```{r da-summary}
res <- fread(params$da_results)

n_total <- nrow(res)
n_sig   <- sum(res$significant, na.rm = TRUE)
n_up    <- sum(res$direction == "Up",   na.rm = TRUE)
n_down  <- sum(res$direction == "Down", na.rm = TRUE)

summary_tbl <- data.table(
  Metric  = c("Total proteins tested",
              sprintf("Significant (FDR<%.2f, |log2FC|>=%.1f)",
                      params$fdr_cutoff, params$fc_cutoff),
              sprintf("Up in %s",   params$condition_b),
              sprintf("Down in %s", params$condition_b)),
  Count   = c(n_total, n_sig, n_up, n_down)
)

kable(summary_tbl) |>
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)
```

### 4.3 Top significant proteins

```{r top-proteins}
sig_proteins <- res[significant == TRUE][order(adj_pvalue)]

if (nrow(sig_proteins) == 0) {
  cat("No significant proteins found at the specified thresholds.")
} else {
  display_cols <- intersect(
    c("gene_names", "protein_id", "log2FC", "pvalue", "adj_pvalue", "direction"),
    colnames(sig_proteins)
  )
  top_n <- head(sig_proteins[, ..display_cols], 30)
  top_n[, log2FC    := round(log2FC, 3)]
  top_n[, pvalue    := signif(pvalue, 3)]
  top_n[, adj_pvalue := signif(adj_pvalue, 3)]

  kable(top_n,
        caption = sprintf("Top %d significant proteins (sorted by adjusted p-value)",
                          nrow(top_n))) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) |>
    row_spec(which(top_n$direction == "Up"),   background = "#fde8e8") |>
    row_spec(which(top_n$direction == "Down"), background = "#e8effe")
}
```

---

## 5. Reproducibility

This analysis was performed with the following software:

```{r session-info}
si <- sessionInfo()
cat(sprintf("R version: %s\n", si$R.version$version.string))
pkgs <- c("limma", "DEqMS", "data.table", "ggplot2")
for (pkg in pkgs) {
  v <- tryCatch(as.character(packageVersion(pkg)), error = function(e) "not installed")
  cat(sprintf("  %s: %s\n", pkg, v))
}
```

> **To reproduce this analysis**, run:
> ```bash
> nextflow run main.nf -profile docker --data_url <URL> --sample_sheet <FILE> \\
>     --condition_a <A> --condition_b <B>
> ```

'

# Write temp Rmd
tmp_rmd <- tempfile(fileext = ".Rmd")
writeLines(rmd_template, tmp_rmd)

cat("[REPORT] Rendering HTML report\n")

rmarkdown::render(
  input       = tmp_rmd,
  output_file = basename(opt$output),
  output_dir  = dirname(normalizePath(opt$output, mustWork = FALSE)),
  params      = list(
    qc_summary   = normalizePath(opt$qc_summary, mustWork = TRUE),
    boxplot_path = normalizePath(opt$boxplot,    mustWork = TRUE),
    da_results   = normalizePath(opt$da_results, mustWork = TRUE),
    volcano_path = normalizePath(opt$volcano,    mustWork = TRUE),
    sample_sheet = normalizePath(opt$sample_sheet, mustWork = TRUE),
    condition_a  = opt$condition_a,
    condition_b  = opt$condition_b,
    fdr_cutoff   = opt$fdr_cutoff,
    fc_cutoff    = opt$fc_cutoff
  ),
  quiet = FALSE
)

cat(sprintf("[REPORT] Wrote report: %s\n", opt$output))
cat("[REPORT] Done.\n")
