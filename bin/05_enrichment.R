#!/usr/bin/env Rscript
# =============================================================================
# 07_enrichment.R
# GO and KEGG pathway enrichment analysis using clusterProfiler
#
# Steps:
#   1. Read DA results (combined multi-comparison TSV)
#   2. For each comparison, extract significant proteins
#   3. Map UniProt IDs to Entrez gene IDs via org.Hs.eg.db
#   4. Run GO enrichment (Biological Process) and KEGG enrichment
#   5. Write enrichment TSVs and dot plot PDFs
#
# Notes:
#   - Only human proteins will map successfully; yeast/non-human skipped
#   - If too few proteins map, enrichment is skipped with a warning
#   - Requires org.Hs.eg.db for human annotation
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    da_results  = NULL,
    output_dir  = ".",
    fdr_cutoff  = 0.05,
    min_genes   = 3      # min mapped genes to attempt enrichment
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--da-results" = { defaults$da_results <- args[i+1]; i <- i+2 },
      "--output-dir" = { defaults$output_dir <- args[i+1]; i <- i+2 },
      "--fdr-cutoff" = { defaults$fdr_cutoff <- as.numeric(args[i+1]); i <- i+2 },
      "--min-genes"  = { defaults$min_genes  <- as.integer(args[i+1]); i <- i+2 },
      { stop(sprintf("Unknown argument: %s", args[i])) }
    )
  }

  if (is.null(defaults$da_results)) stop("--da-results is required")
  defaults
}

opt <- parse_args(args)
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("[ENRICH] Reading DA results:", opt$da_results, "\n")
res <- fread(opt$da_results, sep = "\t", header = TRUE, data.table = TRUE)

# Check required columns
required_cols <- c("protein_id", "gene_names", "significant",
                   "comparison", "log2FC", "adj_pvalue")
missing <- setdiff(required_cols, colnames(res))
if (length(missing) > 0) {
  stop(sprintf("Missing columns in DA results: %s", paste(missing, collapse = ", ")))
}

comparisons <- unique(res$comparison)
cat(sprintf("[ENRICH] Comparisons found: %s\n", paste(comparisons, collapse = ", ")))

# -----------------------------------------------------------------------------
# Helper: map gene names or UniProt IDs to Entrez IDs
# -----------------------------------------------------------------------------
map_to_entrez <- function(gene_names_vec, protein_ids_vec) {

  # Clean gene names — take first gene name if semicolon-separated
  clean_genes <- gsub(";.*", "", gene_names_vec)
  clean_genes <- trimws(clean_genes)
  clean_genes <- clean_genes[clean_genes != "" & !is.na(clean_genes)]
  clean_genes <- unique(clean_genes)

  if (length(clean_genes) == 0) {
    cat("[ENRICH]   No gene names available, trying UniProt IDs\n")
    # Fall back to UniProt IDs
    uniprot_ids <- gsub(";.*", "", protein_ids_vec)
    uniprot_ids <- trimws(uniprot_ids)
    uniprot_ids <- unique(uniprot_ids[uniprot_ids != "" & !is.na(uniprot_ids)])

    mapped <- tryCatch(
      bitr(uniprot_ids, fromType = "UNIPROT",
           toType = "ENTREZID", OrgDb = org.Hs.eg.db),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(
        bitr(uniprot_ids, fromType = "UNIPROT",
             toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      )
    )
    return(mapped$ENTREZID)
  }

  mapped <- tryCatch(
    suppressMessages(suppressWarnings(
      bitr(clean_genes, fromType = "SYMBOL",
           toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    )),
    error = function(e) NULL
  )

  if (is.null(mapped) || nrow(mapped) == 0) return(character(0))
  unique(mapped$ENTREZID)
}

# -----------------------------------------------------------------------------
# Helper: run GO and KEGG enrichment for one gene set
# -----------------------------------------------------------------------------
run_enrichment <- function(entrez_ids, universe_entrez, label, output_dir,
                           fdr_cutoff, min_genes) {

  if (length(entrez_ids) < min_genes) {
    cat(sprintf("[ENRICH]   Skipping %s: only %d genes mapped (need >=%d)\n",
                label, length(entrez_ids), min_genes))
    return(NULL)
  }

  results <- list()

  # --- GO Biological Process ---
  cat(sprintf("[ENRICH]   Running GO BP for %s (%d genes)\n",
              label, length(entrez_ids)))

  go_result <- tryCatch(
    suppressMessages(enrichGO(
      gene          = entrez_ids,
      universe      = universe_entrez,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = fdr_cutoff,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )),
    error = function(e) {
      cat(sprintf("[ENRICH]   GO BP failed: %s\n", e$message))
      NULL
    }
  )

  if (!is.null(go_result) && nrow(go_result@result) > 0) {
    go_df <- as.data.table(go_result@result)
    go_df[, comparison := label]
    go_df[, analysis   := "GO_BP"]
    results[["go"]] <- go_df

    # Dot plot — top 20 terms
    go_plot_path <- file.path(output_dir,
                              sprintf("go_dotplot_%s.pdf", label))
    n_show <- min(20, nrow(go_result@result[go_result@result$p.adjust < fdr_cutoff, ]))
    if (n_show > 0) {
      p <- dotplot(go_result, showCategory = n_show, font.size = 9) +
        labs(title = sprintf("GO Biological Process: %s", gsub("_", " ", label))) +
        theme(plot.title = element_text(size = 11))
      pdf(go_plot_path, width = 10, height = max(5, n_show * 0.4 + 2))
      print(p)
      dev.off()
      cat(sprintf("[ENRICH]   Wrote GO plot: %s\n", go_plot_path))
    }
  } else {
    cat(sprintf("[ENRICH]   No significant GO terms for %s\n", label))
  }

  # --- KEGG ---
  cat(sprintf("[ENRICH]   Running KEGG for %s\n", label))

  kegg_result <- tryCatch(
    suppressMessages(enrichKEGG(
      gene          = entrez_ids,
      universe      = universe_entrez,
      organism      = "hsa",
      pAdjustMethod = "BH",
      pvalueCutoff  = fdr_cutoff,
      qvalueCutoff  = 0.2
    )),
    error = function(e) {
      cat(sprintf("[ENRICH]   KEGG failed: %s\n", e$message))
      NULL
    }
  )

  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    kegg_df <- as.data.table(kegg_result@result)
    kegg_df[, comparison := label]
    kegg_df[, analysis   := "KEGG"]
    results[["kegg"]] <- kegg_df

    n_show <- min(20, nrow(kegg_result@result[kegg_result@result$p.adjust < fdr_cutoff, ]))
    if (n_show > 0) {
      kegg_plot_path <- file.path(output_dir,
                                  sprintf("kegg_dotplot_%s.pdf", label))
      p <- dotplot(kegg_result, showCategory = n_show, font.size = 9) +
        labs(title = sprintf("KEGG Pathways: %s", gsub("_", " ", label))) +
        theme(plot.title = element_text(size = 11))
      pdf(kegg_plot_path, width = 10, height = max(5, n_show * 0.4 + 2))
      print(p)
      dev.off()
      cat(sprintf("[ENRICH]   Wrote KEGG plot: %s\n", kegg_plot_path))
    }
  } else {
    cat(sprintf("[ENRICH]   No significant KEGG pathways for %s\n", label))
  }

  if (length(results) == 0) return(NULL)
  rbindlist(results, use.names = TRUE, fill = TRUE)
}

# -----------------------------------------------------------------------------
# Build universe: all tested proteins mapped to Entrez
# -----------------------------------------------------------------------------
cat("[ENRICH] Building background gene universe\n")
universe_entrez <- map_to_entrez(res$gene_names, res$protein_id)
cat(sprintf("[ENRICH] Universe: %d Entrez IDs mapped\n", length(universe_entrez)))

if (length(universe_entrez) == 0) {
  cat("[ENRICH] WARNING: No genes could be mapped to Entrez IDs.\n")
  cat("[ENRICH] This is expected for non-human datasets (e.g. yeast).\n")
  cat("[ENRICH] Writing empty enrichment file and exiting.\n")

  empty <- data.table(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    comparison = character(), analysis = character()
  )
  fwrite(empty, file.path(opt$output_dir, "enrichment_results.tsv"),
         sep = "\t", quote = FALSE)
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Run enrichment per comparison
# -----------------------------------------------------------------------------
all_enrichment <- list()

for (comp in comparisons) {
  cat(sprintf("[ENRICH] Processing comparison: %s\n", comp))

  sig_proteins <- res[comparison == comp & significant == TRUE]
  cat(sprintf("[ENRICH]   Significant proteins: %d\n", nrow(sig_proteins)))

  if (nrow(sig_proteins) == 0) {
    cat(sprintf("[ENRICH]   Skipping %s: no significant proteins\n", comp))
    next
  }

  entrez_ids <- map_to_entrez(sig_proteins$gene_names, sig_proteins$protein_id)
  cat(sprintf("[ENRICH]   Mapped to Entrez: %d\n", length(entrez_ids)))

  enrich_res <- run_enrichment(
    entrez_ids    = entrez_ids,
    universe_entrez = universe_entrez,
    label         = comp,
    output_dir    = opt$output_dir,
    fdr_cutoff    = opt$fdr_cutoff,
    min_genes     = opt$min_genes
  )

  if (!is.null(enrich_res)) {
    all_enrichment[[comp]] <- enrich_res
  }
}

# -----------------------------------------------------------------------------
# Write combined enrichment results
# -----------------------------------------------------------------------------
out_path <- file.path(opt$output_dir, "enrichment_results.tsv")

if (length(all_enrichment) > 0) {
  combined <- rbindlist(all_enrichment, use.names = TRUE, fill = TRUE)
  fwrite(combined, out_path, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("[ENRICH] Wrote enrichment results: %s (%d rows)\n",
              out_path, nrow(combined)))
} else {
  cat("[ENRICH] No enrichment results to write (no significant terms found).\n")
  empty <- data.table(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    comparison = character(), analysis = character()
  )
  fwrite(empty, out_path, sep = "\t", quote = FALSE)
}

cat("[ENRICH] Done.\n")
