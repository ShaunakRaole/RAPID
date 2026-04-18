#!/usr/bin/env Rscript
# =============================================================================
# 02_normalize.R
# Log2 transformation and median centering normalization
#
# Steps:
#   1. Read filtered matrix from QC step
#   2. Log2 transform intensities
#   3. Median centering (subtract per-sample median from each sample)
#   4. Write normalized matrix
#   5. Write before/after boxplot PDF
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    input        = NULL,
    sample_sheet = NULL,
    output       = "normalized_matrix.tsv",
    boxplot      = "normalization_boxplot.pdf"
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--input"        = { defaults$input        <- args[i+1]; i <- i+2 },
      "--sample-sheet" = { defaults$sample_sheet <- args[i+1]; i <- i+2 },
      "--output"       = { defaults$output        <- args[i+1]; i <- i+2 },
      "--boxplot"      = { defaults$boxplot        <- args[i+1]; i <- i+2 },
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

cat("[NORM] Reading filtered matrix:", opt$input, "\n")
mat <- fread(opt$input, sep = "\t", header = TRUE, data.table = TRUE)
ss  <- fread(opt$sample_sheet, sep = ",", header = TRUE, data.table = TRUE)

intensity_cols <- ss$sample

# Verify all sample columns are present
missing <- setdiff(intensity_cols, colnames(mat))
if (length(missing) > 0) {
  stop(sprintf("Sample columns not found in matrix: %s", paste(missing, collapse = ", ")))
}

# -----------------------------------------------------------------------------
# Helper: melt matrix to long form for ggplot
# -----------------------------------------------------------------------------
melt_intensities <- function(mat, intensity_cols, ss, label) {
  long <- melt(mat,
    id.vars      = c("protein_id", "gene_names"),
    measure.vars = intensity_cols,
    variable.name = "sample",
    value.name   = "intensity"
  )
  long <- long[!is.na(intensity)]
  long <- merge(long, ss[, .(sample, condition)], by = "sample")
  long[, stage := label]
  long
}

# -----------------------------------------------------------------------------
# Step 1: Log2 transform
# -----------------------------------------------------------------------------
cat("[NORM] Applying log2 transformation\n")

mat_log2 <- copy(mat)
for (col in intensity_cols) {
  mat_log2[, (col) := log2(get(col))]
}

# Collect pre-normalization data for plot
long_before <- melt_intensities(mat_log2, intensity_cols, ss, "Before normalization")

# -----------------------------------------------------------------------------
# Step 2: Median centering
# -----------------------------------------------------------------------------
cat("[NORM] Applying median centering\n")

# Calculate global median across all samples (excluding NA)
all_vals <- unlist(mat_log2[, ..intensity_cols], use.names = FALSE)
global_median <- median(all_vals, na.rm = TRUE)
cat(sprintf("[NORM] Global median: %.4f\n", global_median))

mat_norm <- copy(mat_log2)
for (col in intensity_cols) {
  sample_median <- median(mat_norm[[col]], na.rm = TRUE)
  shift <- global_median - sample_median
  cat(sprintf("[NORM]   %s: median=%.4f, shift=%.4f\n", col, sample_median, shift))
  mat_norm[, (col) := get(col) + shift]
}

# Collect post-normalization data for plot
long_after <- melt_intensities(mat_norm, intensity_cols, ss, "After normalization")

# -----------------------------------------------------------------------------
# Step 3: Boxplot
# -----------------------------------------------------------------------------
cat("[NORM] Generating normalization boxplot\n")

long_both <- rbindlist(list(long_before, long_after))
long_both[, stage := factor(stage, levels = c("Before normalization", "After normalization"))]

# Condition colours — use up to 8 colours
cond_colors <- c("#4E9AF1", "#E05C5C", "#50C878", "#FFA500",
                 "#9B59B6", "#1ABC9C", "#E67E22", "#2C3E50")
conditions  <- unique(ss$condition)
color_map   <- setNames(cond_colors[seq_along(conditions)], conditions)

p <- ggplot(long_both, aes(x = sample, y = intensity, fill = condition)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~stage, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = color_map) +
  labs(
    title = "Protein intensity distribution before and after normalization",
    x     = "Sample",
    y     = "log2 LFQ intensity",
    fill  = "Condition"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold"),
    legend.position = "right"
  )

pdf(opt$boxplot, width = max(8, length(intensity_cols) * 0.8 + 2), height = 8)
print(p)
dev.off()

cat(sprintf("[NORM] Wrote boxplot: %s\n", opt$boxplot))

# -----------------------------------------------------------------------------
# Step 4: Write normalized matrix
# -----------------------------------------------------------------------------
fwrite(mat_norm, opt$output, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[NORM] Wrote normalized matrix: %s\n", opt$output))
cat("[NORM] Done.\n")
