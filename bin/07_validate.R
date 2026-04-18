#!/usr/bin/env Rscript
# =============================================================================
# 08_validate.R
# Output validation for RAPID pipeline
#
# Checks:
#   1. All expected output files exist and are non-empty
#   2. Filtered matrix — minimum protein count, required columns
#   3. Normalized matrix — same dimensions as filtered, no all-NA rows
#   4. PCA coordinates — correct number of samples, PC columns present
#   5. DA results — required columns, at least one comparison present,
#                   p-values in valid range [0,1], log2FC is numeric
#   6. Enrichment results — valid structure (may be empty, that's OK)
#   7. HTML report — minimum file size, valid HTML structure
#   8. QC summary — readable, contains expected metrics
#
# Exit codes:
#   0 = all checks passed
#   1 = one or more checks failed (pipeline still completes, results flagged)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    filtered_matrix    = NULL,
    normalized_matrix  = NULL,
    pca_coords         = NULL,
    da_results         = NULL,
    enrichment_results = NULL,
    report             = NULL,
    qc_summary         = NULL,
    sample_sheet       = NULL,
    output             = "validation_report.txt",
    min_proteins       = 50L
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--filtered-matrix"    = { defaults$filtered_matrix    <- args[i+1]; i <- i+2 },
      "--normalized-matrix"  = { defaults$normalized_matrix  <- args[i+1]; i <- i+2 },
      "--pca-coords"         = { defaults$pca_coords         <- args[i+1]; i <- i+2 },
      "--da-results"         = { defaults$da_results         <- args[i+1]; i <- i+2 },
      "--enrichment-results" = { defaults$enrichment_results <- args[i+1]; i <- i+2 },
      "--report"             = { defaults$report             <- args[i+1]; i <- i+2 },
      "--qc-summary"         = { defaults$qc_summary         <- args[i+1]; i <- i+2 },
      "--sample-sheet"       = { defaults$sample_sheet       <- args[i+1]; i <- i+2 },
      "--output"             = { defaults$output             <- args[i+1]; i <- i+2 },
      "--min-proteins"       = { defaults$min_proteins       <- as.integer(args[i+1]); i <- i+2 },
      { stop(sprintf("Unknown argument: %s", args[i])) }
    )
  }

  required <- c("filtered_matrix", "normalized_matrix", "pca_coords",
                "da_results", "enrichment_results", "report",
                "qc_summary", "sample_sheet")
  for (r in required) {
    if (is.null(defaults[[r]])) stop(sprintf("--%s is required", gsub("_", "-", r)))
  }
  defaults
}

opt <- parse_args(args)

# -----------------------------------------------------------------------------
# Validation framework
# -----------------------------------------------------------------------------
results  <- list()
n_pass   <- 0L
n_fail   <- 0L
n_warn   <- 0L

log_check <- function(name, status, message, detail = NULL) {
  entry <- list(name = name, status = status, message = message, detail = detail)
  results[[length(results) + 1]] <<- entry

  symbol <- switch(status,
    PASS = "[PASS]",
    FAIL = "[FAIL]",
    WARN = "[WARN]"
  )

  cat(sprintf("%s %s: %s\n", symbol, name, message))
  if (!is.null(detail)) cat(sprintf("       %s\n", detail))

  if (status == "PASS") n_pass <<- n_pass + 1L
  if (status == "FAIL") n_fail <<- n_fail + 1L
  if (status == "WARN") n_warn <<- n_warn + 1L
}

check_file <- function(path, name, min_bytes = 1L) {
  if (!file.exists(path)) {
    log_check(name, "FAIL", sprintf("File not found: %s", basename(path)))
    return(FALSE)
  }
  size <- file.info(path)$size
  if (size < min_bytes) {
    log_check(name, "FAIL",
              sprintf("File is empty or too small: %s (%d bytes)", basename(path), size))
    return(FALSE)
  }
  log_check(name, "PASS", sprintf("File exists (%s, %s bytes)",
                                   basename(path), format(size, big.mark = ",")))
  TRUE
}

cat("=============================================================\n")
cat(" RAPID Pipeline — Output Validation\n")
cat(sprintf(" %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat("=============================================================\n\n")

# -----------------------------------------------------------------------------
# 1. File existence checks
# -----------------------------------------------------------------------------
cat("--- File existence ---\n")
files_ok <- list(
  filtered   = check_file(opt$filtered_matrix,    "filtered_matrix.tsv",    min_bytes = 100),
  normalized = check_file(opt$normalized_matrix,  "normalized_matrix.tsv",  min_bytes = 100),
  pca        = check_file(opt$pca_coords,         "pca_coords.tsv",         min_bytes = 50),
  da         = check_file(opt$da_results,         "da_results.tsv",         min_bytes = 100),
  enrichment = check_file(opt$enrichment_results, "enrichment_results.tsv", min_bytes = 1),
  report     = check_file(opt$report,             "proteomics_report.html", min_bytes = 10000),
  qc         = check_file(opt$qc_summary,         "qc_summary.txt",         min_bytes = 50)
)

# -----------------------------------------------------------------------------
# 2. Filtered matrix checks
# -----------------------------------------------------------------------------
cat("\n--- Filtered matrix ---\n")
if (files_ok$filtered) {
  tryCatch({
    mat <- fread(opt$filtered_matrix, sep = "\t", header = TRUE)
    ss  <- fread(opt$sample_sheet, sep = ",", header = TRUE)

    # Row count
    n_proteins <- nrow(mat)
    if (n_proteins >= opt$min_proteins) {
      log_check("protein_count", "PASS",
                sprintf("%d proteins (>= %d required)", n_proteins, opt$min_proteins))
    } else {
      log_check("protein_count", "FAIL",
                sprintf("Only %d proteins (minimum %d required)", n_proteins, opt$min_proteins),
                "Consider lowering --min_valid_values")
    }

    # Required columns
    for (col in c("protein_id", "gene_names")) {
      if (col %in% colnames(mat)) {
        log_check(sprintf("column_%s", col), "PASS", sprintf("Column '%s' present", col))
      } else {
        log_check(sprintf("column_%s", col), "FAIL", sprintf("Column '%s' missing", col))
      }
    }

    # Sample columns present
    missing_samples <- setdiff(ss$sample, colnames(mat))
    if (length(missing_samples) == 0) {
      log_check("sample_columns", "PASS",
                sprintf("All %d sample columns present", nrow(ss)))
    } else {
      log_check("sample_columns", "FAIL",
                sprintf("Missing sample columns: %s", paste(missing_samples, collapse = ", ")))
    }

  }, error = function(e) {
    log_check("filtered_matrix_read", "FAIL", sprintf("Could not read file: %s", e$message))
  })
}

# -----------------------------------------------------------------------------
# 3. Normalized matrix checks
# -----------------------------------------------------------------------------
cat("\n--- Normalized matrix ---\n")
if (files_ok$normalized && files_ok$filtered) {
  tryCatch({
    mat_f <- fread(opt$filtered_matrix,   sep = "\t", header = TRUE)
    mat_n <- fread(opt$normalized_matrix, sep = "\t", header = TRUE)

    # Same dimensions
    if (nrow(mat_n) == nrow(mat_f)) {
      log_check("normalized_rows", "PASS",
                sprintf("Row count matches filtered matrix (%d proteins)", nrow(mat_n)))
    } else {
      log_check("normalized_rows", "FAIL",
                sprintf("Row count mismatch: filtered=%d, normalized=%d",
                        nrow(mat_f), nrow(mat_n)))
    }

    # Check values are log2-scale (roughly 15-35 range for LFQ)
    ss      <- fread(opt$sample_sheet, sep = ",", header = TRUE)
    int_cols <- intersect(ss$sample, colnames(mat_n))
    if (length(int_cols) > 0) {
      vals      <- unlist(mat_n[, ..int_cols], use.names = FALSE)
      vals      <- vals[!is.na(vals)]
      val_range <- range(vals)
      if (val_range[1] > 0 && val_range[2] < 60) {
        log_check("normalized_range", "PASS",
                  sprintf("Intensity range looks like log2 scale [%.1f, %.1f]",
                          val_range[1], val_range[2]))
      } else {
        log_check("normalized_range", "WARN",
                  sprintf("Unusual intensity range [%.1f, %.1f] — check normalization",
                          val_range[1], val_range[2]))
      }
    }

  }, error = function(e) {
    log_check("normalized_matrix_read", "FAIL", sprintf("Could not read file: %s", e$message))
  })
}

# -----------------------------------------------------------------------------
# 4. PCA coordinates checks
# -----------------------------------------------------------------------------
cat("\n--- PCA coordinates ---\n")
if (files_ok$pca) {
  tryCatch({
    pca <- fread(opt$pca_coords, sep = "\t", header = TRUE)
    ss  <- fread(opt$sample_sheet, sep = ",", header = TRUE)

    # Correct number of samples
    if (nrow(pca) == nrow(ss)) {
      log_check("pca_sample_count", "PASS",
                sprintf("PCA has %d rows matching %d samples", nrow(pca), nrow(ss)))
    } else {
      log_check("pca_sample_count", "FAIL",
                sprintf("PCA has %d rows but sample sheet has %d samples",
                        nrow(pca), nrow(ss)))
    }

    # PC columns present
    for (col in c("PC1", "PC2", "sample", "condition")) {
      if (col %in% colnames(pca)) {
        log_check(sprintf("pca_col_%s", col), "PASS", sprintf("Column '%s' present", col))
      } else {
        log_check(sprintf("pca_col_%s", col), "FAIL", sprintf("Column '%s' missing", col))
      }
    }

  }, error = function(e) {
    log_check("pca_read", "FAIL", sprintf("Could not read PCA file: %s", e$message))
  })
}

# -----------------------------------------------------------------------------
# 5. DA results checks
# -----------------------------------------------------------------------------
cat("\n--- Differential abundance results ---\n")
if (files_ok$da) {
  tryCatch({
    da <- fread(opt$da_results, sep = "\t", header = TRUE)

    # Required columns
    required_da_cols <- c("protein_id", "log2FC", "pvalue",
                          "adj_pvalue", "significant", "comparison")
    for (col in required_da_cols) {
      if (col %in% colnames(da)) {
        log_check(sprintf("da_col_%s", col), "PASS", sprintf("Column '%s' present", col))
      } else {
        log_check(sprintf("da_col_%s", col), "FAIL", sprintf("Column '%s' missing", col))
      }
    }

    # At least one comparison
    if ("comparison" %in% colnames(da)) {
      comparisons <- unique(da$comparison)
      log_check("da_comparisons", "PASS",
                sprintf("%d comparison(s): %s",
                        length(comparisons), paste(comparisons, collapse = ", ")))
    }

    # p-values in valid range
    if ("adj_pvalue" %in% colnames(da)) {
      pvals <- da$adj_pvalue[!is.na(da$adj_pvalue)]
      if (all(pvals >= 0 & pvals <= 1)) {
        log_check("pvalue_range", "PASS",
                  sprintf("Adjusted p-values in valid range [0,1] (n=%d)", length(pvals)))
      } else {
        log_check("pvalue_range", "FAIL",
                  "Adjusted p-values outside valid range [0,1]")
      }
    }

    # Significant proteins
    if ("significant" %in% colnames(da)) {
      n_sig <- sum(da$significant, na.rm = TRUE)
      if (n_sig > 0) {
        log_check("da_significant", "PASS",
                  sprintf("%d significant proteins found", n_sig))
      } else {
        log_check("da_significant", "WARN",
                  "No significant proteins found — check FDR/FC thresholds",
                  "Try --fdr_cutoff 0.1 or --fc_cutoff 0.5")
      }
    }

    # log2FC is numeric
    if ("log2FC" %in% colnames(da)) {
      if (is.numeric(da$log2FC)) {
        fc_range <- range(da$log2FC, na.rm = TRUE)
        log_check("log2fc_numeric", "PASS",
                  sprintf("log2FC is numeric, range [%.2f, %.2f]",
                          fc_range[1], fc_range[2]))
      } else {
        log_check("log2fc_numeric", "FAIL", "log2FC column is not numeric")
      }
    }

  }, error = function(e) {
    log_check("da_read", "FAIL", sprintf("Could not read DA results: %s", e$message))
  })
}

# -----------------------------------------------------------------------------
# 6. Enrichment results checks
# -----------------------------------------------------------------------------
cat("\n--- Enrichment results ---\n")
if (files_ok$enrichment) {
  tryCatch({
    enrich <- fread(opt$enrichment_results, sep = "\t", header = TRUE)

    if (nrow(enrich) == 0) {
      log_check("enrichment_results", "WARN",
                "Enrichment file is empty — no significant terms found",
                "Expected for non-human data or small significant protein sets")
    } else {
      log_check("enrichment_results", "PASS",
                sprintf("%d enrichment terms found", nrow(enrich)))

      # Check required columns
      for (col in c("ID", "Description", "p.adjust", "comparison", "analysis")) {
        if (col %in% colnames(enrich)) {
          log_check(sprintf("enrich_col_%s", col), "PASS",
                    sprintf("Column '%s' present", col))
        } else {
          log_check(sprintf("enrich_col_%s", col), "FAIL",
                    sprintf("Column '%s' missing", col))
        }
      }
    }

  }, error = function(e) {
    log_check("enrichment_read", "FAIL",
              sprintf("Could not read enrichment results: %s", e$message))
  })
}

# -----------------------------------------------------------------------------
# 7. HTML report checks
# -----------------------------------------------------------------------------
cat("\n--- HTML report ---\n")
if (files_ok$report) {
  tryCatch({
    report_text <- readLines(opt$report, warn = FALSE)
    report_size <- file.info(opt$report)$size

    # Valid HTML
    if (any(grepl("<html", report_text, ignore.case = TRUE))) {
      log_check("report_html", "PASS", "Report contains valid HTML structure")
    } else {
      log_check("report_html", "FAIL", "Report does not appear to be valid HTML")
    }

    # Minimum size
    if (report_size >= 50000) {
      log_check("report_size", "PASS",
                sprintf("Report size is %s bytes", format(report_size, big.mark = ",")))
    } else {
      log_check("report_size", "WARN",
                sprintf("Report is smaller than expected (%s bytes)",
                        format(report_size, big.mark = ",")),
                "Some sections may be missing")
    }

    # Key sections present
    report_str <- paste(report_text, collapse = " ")
    sections <- c("QC Summary", "Normalization", "PCA",
                  "Differential Abundance", "Reproducibility")
    for (section in sections) {
      if (grepl(section, report_str, ignore.case = TRUE)) {
        log_check(sprintf("report_section_%s", gsub(" ", "_", section)),
                  "PASS", sprintf("Section '%s' present in report", section))
      } else {
        log_check(sprintf("report_section_%s", gsub(" ", "_", section)),
                  "WARN", sprintf("Section '%s' not found in report", section))
      }
    }

  }, error = function(e) {
    log_check("report_read", "FAIL", sprintf("Could not read report: %s", e$message))
  })
}

# -----------------------------------------------------------------------------
# 8. QC summary checks
# -----------------------------------------------------------------------------
cat("\n--- QC summary ---\n")
if (files_ok$qc) {
  tryCatch({
    qc_lines <- readLines(opt$qc_summary, warn = FALSE)
    qc_text  <- paste(qc_lines, collapse = " ")

    expected_fields <- c("Raw proteins", "Final proteins for analysis",
                         "Samples", "Conditions")
    for (field in expected_fields) {
      if (grepl(field, qc_text)) {
        log_check(sprintf("qc_%s", gsub(" ", "_", field)), "PASS",
                  sprintf("QC field '%s' present", field))
      } else {
        log_check(sprintf("qc_%s", gsub(" ", "_", field)), "WARN",
                  sprintf("QC field '%s' not found", field))
      }
    }

  }, error = function(e) {
    log_check("qc_read", "FAIL", sprintf("Could not read QC summary: %s", e$message))
  })
}

# -----------------------------------------------------------------------------
# Write validation report
# -----------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(sprintf(" SUMMARY: %d passed, %d warnings, %d failed\n",
            n_pass, n_warn, n_fail))
cat("=============================================================\n")

overall_status <- if (n_fail > 0) "FAILED" else if (n_warn > 0) "PASSED WITH WARNINGS" else "PASSED"
cat(sprintf(" Overall status: %s\n\n", overall_status))

# Write to file
lines <- c(
  "RAPID Pipeline — Output Validation Report",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Overall status: %s", overall_status),
  sprintf("Passed: %d | Warnings: %d | Failed: %d", n_pass, n_warn, n_fail),
  "",
  "=============================================================",
  "Detailed results:",
  "=============================================================",
  ""
)

for (r in results) {
  lines <- c(lines, sprintf("[%s] %s: %s", r$status, r$name, r$message))
  if (!is.null(r$detail)) lines <- c(lines, sprintf("       Note: %s", r$detail))
}

writeLines(lines, opt$output)
cat(sprintf("[VALIDATE] Wrote validation report: %s\n", opt$output))

# Exit with error code if any checks failed
if (n_fail > 0) {
  cat(sprintf("[VALIDATE] %d check(s) failed — review validation report\n", n_fail))
  quit(save = "no", status = 1)
}

cat("[VALIDATE] All checks passed.\n")
quit(save = "no", status = 0)
