#!/usr/bin/env Rscript
# =============================================================================
# 03_da_analysis.R
# Differential abundance analysis using DEqMS
#
# Steps:
#   1. Read normalized matrix
#   2. Build limma design matrix for condition_a vs condition_b
#   3. Fit limma model + eBayes
#   4. Apply DEqMS correction using peptide count
#   5. Write results TSV and volcano plot PDF
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(DEqMS)
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
    condition_a  = NULL,
    condition_b  = NULL,
    fdr_cutoff   = 0.05,
    fc_cutoff    = 1.0,
    output       = "da_results.tsv",
    volcano      = "volcano_plot.pdf"
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--input"        = { defaults$input        <- args[i+1]; i <- i+2 },
      "--sample-sheet" = { defaults$sample_sheet <- args[i+1]; i <- i+2 },
      "--condition-a"  = { defaults$condition_a   <- args[i+1]; i <- i+2 },
      "--condition-b"  = { defaults$condition_b   <- args[i+1]; i <- i+2 },
      "--fdr-cutoff"   = { defaults$fdr_cutoff    <- as.numeric(args[i+1]); i <- i+2 },
      "--fc-cutoff"    = { defaults$fc_cutoff     <- as.numeric(args[i+1]); i <- i+2 },
      "--output"       = { defaults$output        <- args[i+1]; i <- i+2 },
      "--volcano"      = { defaults$volcano       <- args[i+1]; i <- i+2 },
      {
        stop(sprintf("Unknown argument: %s", args[i]))
      }
    )
  }

  if (is.null(defaults$input))       stop("--input is required")
  if (is.null(defaults$sample_sheet))stop("--sample-sheet is required")
  if (is.null(defaults$condition_a)) stop("--condition-a is required")
  if (is.null(defaults$condition_b)) stop("--condition-b is required")
  defaults
}

opt <- parse_args(args)

cat("[DA] Reading normalized matrix:", opt$input, "\n")
mat <- fread(opt$input, sep = "\t", header = TRUE, data.table = TRUE)
ss  <- fread(opt$sample_sheet, sep = ",", header = TRUE, data.table = TRUE)

# Validate conditions exist in sample sheet
all_conditions <- unique(ss$condition)
for (cond in c(opt$condition_a, opt$condition_b)) {
  if (!cond %in% all_conditions) {
    stop(sprintf(
      "Condition '%s' not found in sample sheet. Available: %s",
      cond, paste(all_conditions, collapse = ", ")
    ))
  }
}

# Subset to only the two conditions being compared
ss_sub <- ss[condition %in% c(opt$condition_a, opt$condition_b)]
intensity_cols <- ss_sub$sample

cat(sprintf("[DA] Comparing: %s vs %s\n", opt$condition_a, opt$condition_b))
cat(sprintf("[DA] Samples in comparison: %d\n", nrow(ss_sub)))

# -----------------------------------------------------------------------------
# Build intensity matrix for limma
# -----------------------------------------------------------------------------
intensity_mat <- as.matrix(mat[, ..intensity_cols])
rownames(intensity_mat) <- mat$protein_id

# -----------------------------------------------------------------------------
# Build DEqMS peptide count vector
# DEqMS requires a count of peptides/PSMs per protein — used as a proxy
# for quantification reliability. We use the number of non-NA values per
# protein across samples (a reasonable fallback when peptide count is not
# directly available in the filtered matrix).
# -----------------------------------------------------------------------------
peptide_counts <- rowSums(!is.na(intensity_mat))
peptide_counts[peptide_counts < 1] <- 1  # guard against zeros

# -----------------------------------------------------------------------------
# Limma design matrix
# -----------------------------------------------------------------------------
condition_factor <- factor(
  ss_sub$condition,
  levels = c(opt$condition_a, opt$condition_b)
)

design <- model.matrix(~ 0 + condition_factor)
colnames(design) <- levels(condition_factor)

cat("[DA] Design matrix:\n")
print(design)

# Contrast: condition_b vs condition_a (positive FC = up in B)
contrast_str <- sprintf("%s-%s", opt$condition_b, opt$condition_a)
contrast_mat <- makeContrasts(
  contrasts = contrast_str,
  levels    = design
)

# -----------------------------------------------------------------------------
# Limma fit
# -----------------------------------------------------------------------------
cat("[DA] Fitting limma model\n")
fit <- lmFit(intensity_mat, design)
fit <- contrasts.fit(fit, contrast_mat)
fit <- eBayes(fit)

# -----------------------------------------------------------------------------
# DEqMS correction
# -----------------------------------------------------------------------------
cat("[DA] Applying DEqMS correction\n")
fit$count <- peptide_counts
fit <- spectraCounteBayes(fit)

# -----------------------------------------------------------------------------
# Extract results
# -----------------------------------------------------------------------------
results <- outputResult(fit, coef_col = 1)

# outputResult returns a data frame; convert to data.table
results_dt <- as.data.table(results, keep.rownames = "protein_id")

# Add gene names from original matrix
gene_map <- mat[, .(protein_id, gene_names)]
results_dt <- merge(results_dt, gene_map, by = "protein_id", all.x = TRUE)

# Rename columns for clarity
setnames(results_dt,
  old = c("logFC", "sca.P.Value", "sca.adj.pval"),
  new = c("log2FC", "pvalue", "adj_pvalue"),
  skip_absent = TRUE
)

# Add significance classification
results_dt[, significant := adj_pvalue < opt$fdr_cutoff & abs(log2FC) >= opt$fc_cutoff]
results_dt[, direction   := fcase(
  significant == TRUE & log2FC >  0, "Up",
  significant == TRUE & log2FC <  0, "Down",
  default = "NS"
)]

# Sort by adjusted p-value
setorder(results_dt, adj_pvalue)

n_sig  <- sum(results_dt$significant, na.rm = TRUE)
n_up   <- sum(results_dt$direction == "Up",   na.rm = TRUE)
n_down <- sum(results_dt$direction == "Down", na.rm = TRUE)

cat(sprintf("[DA] Significant proteins (FDR < %.2f, |log2FC| >= %.1f): %d\n",
            opt$fdr_cutoff, opt$fc_cutoff, n_sig))
cat(sprintf("[DA]   Up   (%s): %d\n", opt$condition_b, n_up))
cat(sprintf("[DA]   Down (%s): %d\n", opt$condition_a, n_down))

# -----------------------------------------------------------------------------
# Volcano plot
# -----------------------------------------------------------------------------
cat("[DA] Generating volcano plot\n")

# Label top 15 significant proteins by adjusted p-value
label_dt <- results_dt[significant == TRUE][order(adj_pvalue)][seq_len(min(15, .N))]
label_dt[, label := ifelse(gene_names != "" & !is.na(gene_names),
                           gene_names, protein_id)]

plot_dt <- copy(results_dt)
plot_dt[, neg_log10_p := -log10(pmax(adj_pvalue, 1e-300))]

# Cap extreme values for display
fc_lim <- max(abs(plot_dt$log2FC), na.rm = TRUE) + 0.5
y_lim  <- max(plot_dt$neg_log10_p, na.rm = TRUE) + 1

vp <- ggplot(plot_dt, aes(x = log2FC, y = neg_log10_p, colour = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(opt$fdr_cutoff),
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_vline(xintercept = c(-opt$fc_cutoff, opt$fc_cutoff),
             linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_text_repel(
    data        = label_dt,
    aes(label   = label),
    size        = 3,
    max.overlaps = 20,
    colour      = "black"
  ) +
  scale_colour_manual(
    values = c("Up" = "#E05C5C", "Down" = "#4E9AF1", "NS" = "grey70"),
    labels = c(
      "Up"   = sprintf("Up in %s (%d)", opt$condition_b, n_up),
      "Down" = sprintf("Down in %s (%d)", opt$condition_b, n_down),
      "NS"   = sprintf("Not significant (%d)", sum(plot_dt$direction == "NS", na.rm=TRUE))
    )
  ) +
  labs(
    title    = sprintf("Differential abundance: %s vs %s", opt$condition_b, opt$condition_a),
    subtitle = sprintf("FDR < %.2f, |log2FC| >= %.1f  |  %d significant proteins",
                       opt$fdr_cutoff, opt$fc_cutoff, n_sig),
    x        = sprintf("log2 fold change (%s / %s)", opt$condition_b, opt$condition_a),
    y        = "-log10 adjusted p-value",
    colour   = NULL
  ) +
  xlim(-fc_lim, fc_lim) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

pdf(opt$volcano, width = 8, height = 7)
print(vp)
dev.off()

cat(sprintf("[DA] Wrote volcano plot: %s\n", opt$volcano))

# -----------------------------------------------------------------------------
# Write results
# -----------------------------------------------------------------------------
fwrite(results_dt, opt$output, sep = "\t", quote = FALSE, na = "NA")
cat(sprintf("[DA] Wrote results: %s\n", opt$output))
cat("[DA] Done.\n")
