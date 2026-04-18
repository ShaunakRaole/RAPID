#!/usr/bin/env Rscript
# =============================================================================
# 01_qc_filter.R
# QC filtering of MaxQuant proteinGroups.txt
#
# Steps:
#   1. Read proteinGroups.txt
#   2. Remove reverse hits, contaminants, only-identified-by-site proteins
#   3. Extract LFQ intensity columns matching the sample sheet
#   4. Replace zero intensities with NA
#   5. Remove proteins with fewer than --min-valid valid values across all samples
#   6. Write filtered matrix and QC summary
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
    input        = NULL,
    sample_sheet = NULL,
    min_valid    = 2L,
    output       = "filtered_matrix.tsv",
    qc_summary   = "qc_summary.txt"
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--input"        = { defaults$input        <- args[i+1]; i <- i+2 },
      "--sample-sheet" = { defaults$sample_sheet <- args[i+1]; i <- i+2 },
      "--min-valid"    = { defaults$min_valid     <- as.integer(args[i+1]); i <- i+2 },
      "--output"       = { defaults$output        <- args[i+1]; i <- i+2 },
      "--qc-summary"   = { defaults$qc_summary    <- args[i+1]; i <- i+2 },
      {
        stop(sprintf("Unknown argument: %s", args[i]))
      }
    )
  }

  if (is.null(defaults$input))        stop("--input is required")
  if (is.null(defaults$sample_sheet)) stop("--sample-sheet is required")
  defaults
}

opt <- parse_args(args)

cat("[QC] Reading proteinGroups:", opt$input, "\n")
cat("[QC] Reading sample sheet: ", opt$sample_sheet, "\n")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
pg <- fread(opt$input, sep = "\t", header = TRUE, data.table = TRUE)
ss <- fread(opt$sample_sheet, sep = ",", header = TRUE, data.table = TRUE)

# Validate sample sheet columns
required_cols <- c("sample", "condition", "file_name")
missing <- setdiff(required_cols, colnames(ss))
if (length(missing) > 0) {
  stop(sprintf("Sample sheet missing columns: %s", paste(missing, collapse = ", ")))
}

n_raw <- nrow(pg)
cat(sprintf("[QC] Raw proteins: %d\n", n_raw))

# -----------------------------------------------------------------------------
# Remove standard MaxQuant contaminants / artifacts
# -----------------------------------------------------------------------------

# Column names vary slightly across MaxQuant versions — handle both
reverse_col    <- intersect(c("Reverse", "reverse"), colnames(pg))[1]
contam_col     <- intersect(c("Potential contaminant", "Potential.contaminant",
                               "potential contaminant"), colnames(pg))[1]
site_only_col  <- intersect(c("Only identified by site", "Only.identified.by.site",
                               "only identified by site"), colnames(pg))[1]

n_before <- nrow(pg)

if (!is.na(reverse_col)) {
  pg <- pg[is.na(get(reverse_col)) | get(reverse_col) != "+"]
  cat(sprintf("[QC] After removing reverse hits: %d proteins\n", nrow(pg)))
}

if (!is.na(contam_col)) {
  pg <- pg[is.na(get(contam_col)) | get(contam_col) != "+"]
  cat(sprintf("[QC] After removing contaminants: %d proteins\n", nrow(pg)))
}

if (!is.na(site_only_col)) {
  pg <- pg[is.na(get(site_only_col)) | get(site_only_col) != "+"]
  cat(sprintf("[QC] After removing site-only: %d proteins\n", nrow(pg)))
}

n_after_flags <- nrow(pg)

# -----------------------------------------------------------------------------
# Extract LFQ intensity columns
# -----------------------------------------------------------------------------
# MaxQuant names LFQ columns as "LFQ intensity <sample_file_name>"
lfq_cols <- paste0("LFQ intensity ", ss$file_name)

# Check all expected columns are present
missing_lfq <- setdiff(lfq_cols, colnames(pg))
if (length(missing_lfq) > 0) {
  cat("[QC] WARNING: The following expected LFQ columns were not found:\n")
  cat(paste(" ", missing_lfq, collapse = "\n"), "\n")
  cat("[QC] Available LFQ columns:\n")
  cat(paste(" ", grep("^LFQ intensity", colnames(pg), value = TRUE), collapse = "\n"), "\n")
  stop("LFQ column mismatch between sample sheet and proteinGroups.txt. Check file_name column in sample sheet.")
}

# Protein identifier column
id_col <- intersect(c("Protein IDs", "Protein.IDs", "protein IDs"), colnames(pg))[1]
if (is.na(id_col)) stop("Cannot find 'Protein IDs' column in proteinGroups.txt")

# Gene names column (optional but useful)
gene_col <- intersect(c("Gene names", "Gene.names", "gene names"), colnames(pg))[1]

# Build intensity matrix
if (!is.na(gene_col)) {
  mat <- pg[, c(id_col, gene_col, lfq_cols), with = FALSE]
  setnames(mat, c("protein_id", "gene_names", ss$sample))
} else {
  mat <- pg[, c(id_col, lfq_cols), with = FALSE]
  setnames(mat, c("protein_id", ss$sample))
  mat[, gene_names := protein_id]
}

# -----------------------------------------------------------------------------
# Replace zeros with NA (zero = not detected in MaxQuant LFQ)
# -----------------------------------------------------------------------------
intensity_cols <- ss$sample
for (col in intensity_cols) {
  mat[get(col) == 0, (col) := NA_real_]
}

# -----------------------------------------------------------------------------
# Filter: remove proteins below minimum valid value threshold
# -----------------------------------------------------------------------------
valid_counts <- rowSums(!is.na(mat[, ..intensity_cols]))
mat <- mat[valid_counts >= opt$min_valid]

n_final <- nrow(mat)
cat(sprintf("[QC] After min-valid filter (>= %d): %d proteins\n", opt$min_valid, n_final))

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
fwrite(mat, opt$output, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[QC] Wrote filtered matrix: %s\n", opt$output))

# QC summary
summary_lines <- c(
  sprintf("Raw proteins:                  %d", n_raw),
  sprintf("After flag removal:            %d", n_after_flags),
  sprintf("Removed (flags):               %d", n_raw - n_after_flags),
  sprintf("After min-valid filter (>=%d): %d", opt$min_valid, n_final),
  sprintf("Removed (low coverage):        %d", n_after_flags - n_final),
  sprintf("Final proteins for analysis:   %d", n_final),
  sprintf("Samples:                       %d", length(intensity_cols)),
  sprintf("Conditions:                    %s", paste(unique(ss$condition), collapse = ", "))
)

writeLines(summary_lines, opt$qc_summary)
cat("[QC] Done.\n")
