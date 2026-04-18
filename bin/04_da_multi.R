#!/usr/bin/env Rscript
# =============================================================================
# 06_da_multi.R
# Differential abundance analysis for all pairwise condition comparisons
#
# Steps:
#   1. Read normalized matrix and sample sheet
#   2. Detect all unique conditions
#   3. For each pairwise combination, run limma + DEqMS
#   4. Write combined results TSV and one volcano PDF per comparison
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
    input       = NULL,
    sample_sheet = NULL,
    fdr_cutoff  = 0.05,
    fc_cutoff   = 1.0,
    output_dir  = "."
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--input"        = { defaults$input        <- args[i+1]; i <- i+2 },
      "--sample-sheet" = { defaults$sample_sheet <- args[i+1]; i <- i+2 },
      "--fdr-cutoff"   = { defaults$fdr_cutoff   <- as.numeric(args[i+1]); i <- i+2 },
      "--fc-cutoff"    = { defaults$fc_cutoff    <- as.numeric(args[i+1]); i <- i+2 },
      "--output-dir"   = { defaults$output_dir   <- args[i+1]; i <- i+2 },
      { stop(sprintf("Unknown argument: %s", args[i])) }
    )
  }

  if (is.null(defaults$input))        stop("--input is required")
  if (is.null(defaults$sample_sheet)) stop("--sample-sheet is required")
  defaults
}

opt <- parse_args(args)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("[DA_MULTI] Reading normalized matrix:", opt$input, "\n")
mat <- fread(opt$input, sep = "\t", header = TRUE, data.table = TRUE)
ss  <- fread(opt$sample_sheet, sep = ",", header = TRUE, data.table = TRUE)

conditions <- unique(ss$condition)
n_cond     <- length(conditions)

if (n_cond < 2) stop("Need at least 2 conditions for pairwise comparisons")

pairs <- combn(conditions, 2, simplify = FALSE)
cat(sprintf("[DA_MULTI] Conditions: %s\n", paste(conditions, collapse = ", ")))
cat(sprintf("[DA_MULTI] Pairwise comparisons: %d\n", length(pairs)))

# -----------------------------------------------------------------------------
# Helper: run one pairwise comparison
# -----------------------------------------------------------------------------
run_comparison <- function(mat, ss, cond_a, cond_b, fdr_cutoff, fc_cutoff) {

  ss_sub        <- ss[condition %in% c(cond_a, cond_b)]
  intensity_cols <- ss_sub$sample

  intensity_mat <- as.matrix(mat[, ..intensity_cols])
  rownames(intensity_mat) <- mat$protein_id

  peptide_counts <- rowSums(!is.na(intensity_mat))
  peptide_counts[peptide_counts < 1] <- 1

  safe_a <- make.names(cond_a)
  safe_b <- make.names(cond_b)

  condition_factor <- factor(make.names(ss_sub$condition),
                             levels = c(safe_a, safe_b))
  design           <- model.matrix(~ 0 + condition_factor)
  colnames(design) <- levels(condition_factor)

  contrast_str <- sprintf("%s-%s", safe_b, safe_a)
  contrast_mat <- makeContrasts(contrasts = contrast_str, levels = design)

  fit        <- lmFit(intensity_mat, design)
  fit        <- contrasts.fit(fit, contrast_mat)
  fit        <- eBayes(fit)
  fit$count  <- peptide_counts
  fit        <- spectraCounteBayes(fit)

  results <- outputResult(fit, coef_col = 1)
  results_dt <- as.data.table(results, keep.rownames = "protein_id")

  gene_map   <- mat[, .(protein_id, gene_names)]
  results_dt <- merge(results_dt, gene_map, by = "protein_id", all.x = TRUE)

  setnames(results_dt,
    old = c("logFC", "sca.P.Value", "sca.adj.pval"),
    new = c("log2FC", "pvalue", "adj_pvalue"),
    skip_absent = TRUE
  )

  results_dt[, comparison  := sprintf("%s_vs_%s", cond_b, cond_a)]
  results_dt[, condition_a := cond_a]
  results_dt[, condition_b := cond_b]
  results_dt[, significant := adj_pvalue < fdr_cutoff & abs(log2FC) >= fc_cutoff]
  results_dt[, direction   := fcase(
    significant == TRUE & log2FC >  0, "Up",
    significant == TRUE & log2FC <  0, "Down",
    default = "NS"
  )]

  setorder(results_dt, adj_pvalue)
  results_dt
}

# -----------------------------------------------------------------------------
# Helper: volcano plot for one comparison
# -----------------------------------------------------------------------------
make_volcano <- function(results_dt, cond_a, cond_b, fdr_cutoff, fc_cutoff) {

  n_up   <- sum(results_dt$direction == "Up",   na.rm = TRUE)
  n_down <- sum(results_dt$direction == "Down", na.rm = TRUE)
  n_sig  <- sum(results_dt$significant, na.rm = TRUE)

  plot_dt <- copy(results_dt)
  plot_dt[, neg_log10_p := -log10(pmax(adj_pvalue, 1e-300))]

  label_dt <- plot_dt[significant == TRUE][order(adj_pvalue)][seq_len(min(15, .N))]
  label_dt[, label := ifelse(!is.na(gene_names) & gene_names != "",
                             gene_names, protein_id)]

  fc_lim <- max(abs(plot_dt$log2FC), na.rm = TRUE) + 0.5

  ggplot(plot_dt, aes(x = log2FC, y = neg_log10_p, colour = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff),
               linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_text_repel(data = label_dt, aes(label = label),
                    size = 3, max.overlaps = 20, colour = "black") +
    scale_colour_manual(
      values = c("Up" = "#E05C5C", "Down" = "#4E9AF1", "NS" = "grey70"),
      labels = c(
        "Up"   = sprintf("Up in %s (%d)", cond_b, n_up),
        "Down" = sprintf("Down in %s (%d)", cond_b, n_down),
        "NS"   = sprintf("Not significant (%d)",
                         sum(plot_dt$direction == "NS", na.rm = TRUE))
      )
    ) +
    labs(
      title    = sprintf("Differential abundance: %s vs %s", cond_b, cond_a),
      subtitle = sprintf("FDR < %.2f, |log2FC| >= %.1f  |  %d significant",
                         fdr_cutoff, fc_cutoff, n_sig),
      x        = sprintf("log2 fold change (%s / %s)", cond_b, cond_a),
      y        = "-log10 adjusted p-value",
      colour   = NULL
    ) +
    xlim(-fc_lim, fc_lim) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

# -----------------------------------------------------------------------------
# Run all pairwise comparisons
# -----------------------------------------------------------------------------
all_results <- list()

for (pair in pairs) {
  cond_a <- pair[1]
  cond_b <- pair[2]
  comp   <- sprintf("%s_vs_%s", cond_b, cond_a)

  cat(sprintf("[DA_MULTI] Running: %s\n", comp))

  res <- tryCatch(
    run_comparison(mat, ss, cond_a, cond_b, opt$fdr_cutoff, opt$fc_cutoff),
    error = function(e) {
      cat(sprintf("[DA_MULTI] WARNING: comparison %s failed: %s\n", comp, e$message))
      NULL
    }
  )

  if (is.null(res)) next

  n_sig <- sum(res$significant, na.rm = TRUE)
  cat(sprintf("[DA_MULTI]   Significant proteins: %d\n", n_sig))

  all_results[[comp]] <- res

  # Write per-comparison volcano
  vp       <- make_volcano(res, cond_a, cond_b, opt$fdr_cutoff, opt$fc_cutoff)
  vol_path <- file.path(opt$output_dir, sprintf("volcano_%s.pdf", comp))
  pdf(vol_path, width = 8, height = 7)
  print(vp)
  dev.off()
  cat(sprintf("[DA_MULTI]   Wrote volcano: %s\n", vol_path))

  # Write per-comparison results
  res_path <- file.path(opt$output_dir, sprintf("da_results_%s.tsv", comp))
  fwrite(res, res_path, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("[DA_MULTI]   Wrote results: %s\n", res_path))
}

# -----------------------------------------------------------------------------
# Write combined results across all comparisons
# -----------------------------------------------------------------------------
if (length(all_results) > 0) {
  combined <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
  combined_path <- file.path(opt$output_dir, "da_results_all_comparisons.tsv")
  fwrite(combined, combined_path, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("[DA_MULTI] Wrote combined results: %s\n", combined_path))
} else {
  stop("[DA_MULTI] All comparisons failed. Check input data.")
}

cat("[DA_MULTI] Done.\n")
