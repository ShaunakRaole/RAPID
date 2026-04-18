#!/usr/bin/env Rscript
# =============================================================================
# 05_pca.R
# PCA of normalized protein intensity matrix
#
# Steps:
#   1. Read normalized matrix
#   2. Impute remaining NAs with row means for PCA
#   3. Run PCA on transposed matrix (samples as rows)
#   4. Plot PC1 vs PC2 coloured by condition
#   5. Write PCA coordinates TSV and plot PDF
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    input        = NULL,
    sample_sheet = NULL,
    output_plot  = "pca_plot.pdf",
    output_data  = "pca_coords.tsv"
  )
  
  i <- 1
  while (i <= length(args)) {
    switch(args[i],
           "--input"        = { defaults$input        <- args[i+1]; i <- i+2 },
           "--sample-sheet" = { defaults$sample_sheet <- args[i+1]; i <- i+2 },
           "--output-plot"  = { defaults$output_plot  <- args[i+1]; i <- i+2 },
           "--output-data"  = { defaults$output_data  <- args[i+1]; i <- i+2 },
           { stop(sprintf("Unknown argument: %s", args[i])) }
    )
  }
  
  if (is.null(defaults$input))        stop("--input is required")
  if (is.null(defaults$sample_sheet)) stop("--sample-sheet is required")
  defaults
}

opt <- parse_args(args)

cat("[PCA] Reading normalized matrix:", opt$input, "\n")
mat <- fread(opt$input, sep = "\t", header = TRUE, data.table = TRUE)
ss  <- fread(opt$sample_sheet, sep = ",", header = TRUE, data.table = TRUE)

intensity_cols <- ss$sample

# Verify columns present
missing <- setdiff(intensity_cols, colnames(mat))
if (length(missing) > 0) {
  stop(sprintf("Sample columns not found in matrix: %s", paste(missing, collapse = ", ")))
}

# -----------------------------------------------------------------------------
# Extract intensity matrix and impute NAs with row means
# PCA requires complete data — row-mean imputation is standard for proteomics
# -----------------------------------------------------------------------------
int_mat <- as.matrix(mat[, ..intensity_cols])
rownames(int_mat) <- mat$protein_id

# Row-mean imputation
row_means <- rowMeans(int_mat, na.rm = TRUE)
for (j in seq_len(ncol(int_mat))) {
  na_idx <- is.na(int_mat[, j])
  int_mat[na_idx, j] <- row_means[na_idx]
}

# Remove any rows that are still all NA after imputation (all missing)
valid_rows <- rowSums(is.finite(int_mat)) == ncol(int_mat)
int_mat <- int_mat[valid_rows, ]
cat(sprintf("[PCA] Proteins used for PCA: %d\n", nrow(int_mat)))

# -----------------------------------------------------------------------------
# Run PCA — samples as rows, proteins as columns
# -----------------------------------------------------------------------------
pca_result <- prcomp(t(int_mat), scale. = TRUE, center = TRUE)

var_explained <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

pca_df <- data.table(
  sample = rownames(pca_result$x),
  PC1    = pca_result$x[, 1],
  PC2    = pca_result$x[, 2]
)

# Merge with condition info
pca_df <- merge(pca_df, ss[, .(sample, condition)], by = "sample")

cat(sprintf("[PCA] PC1: %.1f%% variance, PC2: %.1f%% variance\n",
            var_explained[1], var_explained[2]))

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
cond_colors <- c("#4E9AF1", "#E05C5C", "#50C878", "#FFA500",
                 "#9B59B6", "#1ABC9C", "#E67E22", "#2C3E50")
conditions  <- unique(pca_df$condition)
color_map   <- setNames(cond_colors[seq_along(conditions)], conditions)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition, label = sample)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_colour_manual(values = color_map) +
  labs(
    title    = "PCA of normalized protein intensities",
    subtitle = sprintf("Each point is one sample  |  %d proteins", nrow(int_mat)),
    x        = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
    y        = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
    colour   = "Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right")

png(opt$output_plot, width = 800, height = 650, res = 120)
print(p)
dev.off()
cat(sprintf("[PCA] Wrote PCA plot: %s\n", opt$output_plot))

# -----------------------------------------------------------------------------
# Write coordinates
# -----------------------------------------------------------------------------
fwrite(pca_df, opt$output_data, sep = "\t", quote = FALSE)
cat(sprintf("[PCA] Wrote PCA coordinates: %s\n", opt$output_data))
cat("[PCA] Done.\n")