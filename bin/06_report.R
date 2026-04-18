#!/usr/bin/env Rscript
# =============================================================================
# 04_report.R  (updated)
# Generates a self-contained HTML report summarising the pipeline run
# Now includes: PCA, multi-comparison DA results, pathway enrichment
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
    qc_summary         = NULL,
    boxplot            = NULL,
    pca_plot           = NULL,
    da_results         = NULL,
    enrichment_results = NULL,
    sample_sheet       = NULL,
    condition_a        = NULL,
    condition_b        = NULL,
    fdr_cutoff         = 0.05,
    fc_cutoff          = 1.0,
    output             = "proteomics_report.html"
  )
  
  i <- 1
  while (i <= length(args)) {
    switch(args[i],
           "--qc-summary"         = { defaults$qc_summary         <- args[i+1]; i <- i+2 },
           "--boxplot"            = { defaults$boxplot             <- args[i+1]; i <- i+2 },
           "--pca-plot"           = { defaults$pca_plot           <- args[i+1]; i <- i+2 },
           "--da-results"         = { defaults$da_results         <- args[i+1]; i <- i+2 },
           "--enrichment-results" = { defaults$enrichment_results <- args[i+1]; i <- i+2 },
           "--sample-sheet"       = { defaults$sample_sheet       <- args[i+1]; i <- i+2 },
           "--condition-a"        = { defaults$condition_a        <- args[i+1]; i <- i+2 },
           "--condition-b"        = { defaults$condition_b        <- args[i+1]; i <- i+2 },
           "--fdr-cutoff"         = { defaults$fdr_cutoff         <- as.numeric(args[i+1]); i <- i+2 },
           "--fc-cutoff"          = { defaults$fc_cutoff          <- as.numeric(args[i+1]); i <- i+2 },
           "--output"             = { defaults$output             <- args[i+1]; i <- i+2 },
           { stop(sprintf("Unknown argument: %s", args[i])) }
    )
  }
  
  for (req in c("qc_summary", "da_results", "sample_sheet")) {
    if (is.null(defaults[[req]])) {
      stop(sprintf("--%s is required", gsub("_", "-", req)))
    }
  }
  defaults
}

opt <- parse_args(args)

# -----------------------------------------------------------------------------
# Rmd template
# -----------------------------------------------------------------------------
rmd_template <- '---
title: "RAPID — Analysis Report"
date: "`r format(Sys.time(), \'%B %d, %Y\')`"
output:
  html_document:
    theme: flatly
    toc: true
    toc_float: true
    self_contained: true
    code_folding: hide
params:
  qc_summary:         ""
  boxplot_path:       ""
  pca_plot_path:      ""
  da_results:         ""
  enrichment_results: ""
  sample_sheet:       ""
  condition_a:        ""
  condition_b:        ""
  fdr_cutoff:         0.05
  fc_cutoff:          1.0
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,
                      out.width = "100%")
library(data.table)
library(knitr)
library(kableExtra)
library(ggplot2)
```

## Overview

This report summarises a reproducible proteomics differential abundance analysis
produced by **[RAPID](https://github.com/ShaunakRaole/RAPID)** — a Nextflow-based
pipeline implementing:

- **QC filtering**: removal of contaminants, reverse hits, and low-coverage proteins
- **Normalization**: log2 transformation and median centering
- **PCA**: quality control visualization of sample clustering
- **Differential abundance**: DEqMS (all pairwise comparisons)
- **Pathway enrichment**: GO Biological Process and KEGG (clusterProfiler)

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

```{r boxplot, fig.cap="Protein intensity distributions before and after normalization."}
knitr::include_graphics(params$boxplot_path)
```

---

## 4. PCA

PCA provides a global view of sample similarity. Samples from the same condition
should cluster together. Outliers or batch effects become visible here.

```{r pca-plot, fig.cap="PCA of normalized protein intensities. Each point is one sample."}
if (nchar(params$pca_plot_path) > 0 && file.exists(params$pca_plot_path)) {
  knitr::include_graphics(params$pca_plot_path)
} else {
  cat("PCA plot not available.")
}
```

---

## 5. Differential Abundance

All pairwise condition comparisons were tested using DEqMS.

```{r da-load}
res <- fread(params$da_results)
comparisons <- if ("comparison" %in% colnames(res)) unique(res$comparison) else "primary"
```

There are **`r length(comparisons)`** pairwise comparison(s):
**`r paste(comparisons, collapse = ", ")`**.

```{r da-summary-table}
if ("comparison" %in% colnames(res)) {
  summary_tbl <- res[, .(
    Total       = .N,
    Significant = sum(significant, na.rm = TRUE),
    Up          = sum(direction == "Up",   na.rm = TRUE),
    Down        = sum(direction == "Down", na.rm = TRUE)
  ), by = comparison]
  setnames(summary_tbl, "comparison", "Comparison")
  kable(summary_tbl,
        caption = sprintf(
          "DA summary per comparison (FDR < %.2f, |log2FC| >= %.1f)",
          params$fdr_cutoff, params$fc_cutoff)) |>
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)
} else {
  summary_tbl <- data.table(
    Metric = c("Total proteins", "Significant", "Up", "Down"),
    Count  = c(nrow(res),
               sum(res$significant, na.rm = TRUE),
               sum(res$direction == "Up",   na.rm = TRUE),
               sum(res$direction == "Down", na.rm = TRUE))
  )
  kable(summary_tbl) |>
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)
}
```

### Top significant proteins

```{r top-proteins}
sig <- res[significant == TRUE][order(adj_pvalue)]
if (nrow(sig) == 0) {
  cat("No significant proteins found at the specified thresholds.")
} else {
  display_cols <- intersect(
    c("comparison", "gene_names", "protein_id",
      "log2FC", "pvalue", "adj_pvalue", "direction"),
    colnames(sig)
  )
  top_n <- head(sig[, ..display_cols], 30)
  top_n[, log2FC     := round(log2FC, 3)]
  top_n[, pvalue     := signif(pvalue, 3)]
  top_n[, adj_pvalue := signif(adj_pvalue, 3)]
  kable(top_n, caption = sprintf("Top %d significant proteins", nrow(top_n))) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) |>
    row_spec(which(top_n$direction == "Up"),   background = "#fde8e8") |>
    row_spec(which(top_n$direction == "Down"), background = "#e8effe")
}
```

---

## 6. Pathway Enrichment

GO Biological Process and KEGG pathway enrichment was performed on significant
proteins from each comparison using clusterProfiler (human annotation only).

```{r enrichment}
enrich_path <- params$enrichment_results
if (nchar(enrich_path) > 0 && file.exists(enrich_path)) {
  enrich <- fread(enrich_path)
  if (nrow(enrich) == 0) {
    cat("No significant enrichment terms found.\n")
    cat("Note: enrichment requires human proteins mappable to Entrez gene IDs.\n")
    cat("For non-human datasets (e.g. mixed yeast/human), this section may be empty.\n")
  } else {
    enrich_summary <- enrich[p.adjust < 0.05, .(Terms = .N),
                             by = .(comparison, analysis)]
    kable(enrich_summary,
          caption = "Significant enrichment terms per comparison") |>
      kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)

    top_terms <- enrich[p.adjust < 0.05][order(p.adjust)][
      , head(.SD, 10), by = .(comparison, analysis)
    ]
    display_cols <- intersect(
      c("comparison", "analysis", "Description",
        "GeneRatio", "p.adjust", "Count"),
      colnames(top_terms)
    )
    top_terms[, p.adjust := signif(p.adjust, 3)]
    kable(top_terms[, ..display_cols],
          caption = "Top enrichment terms (FDR < 0.05)") |>
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                    full_width = FALSE)
  }
} else {
  cat("Enrichment results not available.\n")
}
```

---

## 7. Reproducibility

```{r session-info}
si <- sessionInfo()
cat(sprintf("R version: %s\n", si$R.version$version.string))
pkgs <- c("limma", "DEqMS", "clusterProfiler", "data.table", "ggplot2")
for (pkg in pkgs) {
  v <- tryCatch(as.character(packageVersion(pkg)),
                error = function(e) "not installed")
  cat(sprintf("  %s: %s\n", pkg, v))
}
```

> **To reproduce this analysis**, run:
> ```bash
> nextflow run main.nf -profile docker \\
>     --data_url    <URL> \\
>     --sample_sheet <FILE> \\
>     --condition_a  <A> \\
>     --condition_b  <B>
> ```

'

tmp_rmd <- tempfile(fileext = ".Rmd")
writeLines(rmd_template, tmp_rmd)

cat("[REPORT] Rendering HTML report\n")
wd <- getwd()

rmarkdown::render(
  input       = tmp_rmd,
  output_file = basename(opt$output),
  output_dir  = wd,
  params      = list(
    qc_summary         = file.path(wd, opt$qc_summary),
    boxplot_path       = file.path(wd, opt$boxplot),
    pca_plot_path      = if (!is.null(opt$pca_plot) && nchar(opt$pca_plot) > 0)
      file.path(wd, opt$pca_plot) else "",
    da_results         = file.path(wd, opt$da_results),
    enrichment_results = if (!is.null(opt$enrichment_results) &&
                             nchar(opt$enrichment_results) > 0)
      file.path(wd, opt$enrichment_results) else "",
    sample_sheet       = file.path(wd, opt$sample_sheet),
    condition_a        = if (!is.null(opt$condition_a)) opt$condition_a else "",
    condition_b        = if (!is.null(opt$condition_b)) opt$condition_b else "",
    fdr_cutoff         = opt$fdr_cutoff,
    fc_cutoff          = opt$fc_cutoff
  ),
  quiet = FALSE
)

cat(sprintf("[REPORT] Wrote report: %s\n", opt$output))
cat("[REPORT] Done.\n")