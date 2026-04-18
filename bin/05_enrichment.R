#!/usr/bin/env Rscript
# =============================================================================
# 07_enrichment.R
# GO and KEGG pathway enrichment — organism-aware
#
# Supported organisms (--organism flag):
#   human  -> org.Hs.eg.db  / KEGG: hsa
#   mouse  -> org.Mm.eg.db  / KEGG: mmu
#   rat    -> org.Rn.eg.db  / KEGG: rno
#   yeast  -> org.Sc.sgd.db / KEGG: sce
#
# Steps:
#   1. Read DA results (combined multi-comparison TSV)
#   2. For each comparison, extract significant proteins
#   3. Map gene names / UniProt IDs to Entrez IDs using the correct OrgDb
#   4. Run GO (Biological Process) and KEGG enrichment
#   5. Write enrichment TSVs and dot plot PDFs
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
})

# -----------------------------------------------------------------------------
# Organism lookup table
# -----------------------------------------------------------------------------
ORGANISM_MAP <- list(
  human = list(
    pkg      = "org.Hs.eg.db",
    kegg     = "hsa",
    id_type  = "SYMBOL",
    fallback = "UNIPROT"
  ),
  mouse = list(
    pkg      = "org.Mm.eg.db",
    kegg     = "mmu",
    id_type  = "SYMBOL",
    fallback = "UNIPROT"
  ),
  rat = list(
    pkg      = "org.Rn.eg.db",
    kegg     = "rno",
    id_type  = "SYMBOL",
    fallback = "UNIPROT"
  ),
  yeast = list(
    pkg      = "org.Sc.sgd.db",
    kegg     = "sce",
    id_type  = "ORF",        # yeast uses ORF names not gene symbols
    fallback = "GENENAME"
  )
)

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  defaults <- list(
    da_results  = NULL,
    output_dir  = ".",
    organism    = "human",
    fdr_cutoff  = 0.05,
    min_genes   = 3
  )

  i <- 1
  while (i <= length(args)) {
    switch(args[i],
      "--da-results" = { defaults$da_results <- args[i+1]; i <- i+2 },
      "--output-dir" = { defaults$output_dir <- args[i+1]; i <- i+2 },
      "--organism"   = { defaults$organism   <- tolower(args[i+1]); i <- i+2 },
      "--fdr-cutoff" = { defaults$fdr_cutoff <- as.numeric(args[i+1]); i <- i+2 },
      "--min-genes"  = { defaults$min_genes  <- as.integer(args[i+1]); i <- i+2 },
      { stop(sprintf("Unknown argument: %s", args[i])) }
    )
  }

  if (is.null(defaults$da_results)) stop("--da-results is required")

  if (!defaults$organism %in% names(ORGANISM_MAP)) {
    stop(sprintf(
      "Unsupported organism '%s'. Choose from: %s",
      defaults$organism,
      paste(names(ORGANISM_MAP), collapse = ", ")
    ))
  }

  defaults
}

opt      <- parse_args(args)
org_info <- ORGANISM_MAP[[opt$organism]]

cat(sprintf("[ENRICH] Organism: %s\n", opt$organism))
cat(sprintf("[ENRICH] OrgDb:    %s\n", org_info$pkg))
cat(sprintf("[ENRICH] KEGG:     %s\n", org_info$kegg))

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load the appropriate OrgDb package dynamically
# -----------------------------------------------------------------------------
if (!requireNamespace(org_info$pkg, quietly = TRUE)) {
  stop(sprintf(
    "Package '%s' is not installed. Add it to your environment.yml and rebuild.\n%s",
    org_info$pkg,
    sprintf("  BiocManager::install('%s')", org_info$pkg)
  ))
}

orgdb <- get(org_info$pkg,
             envir = loadNamespace(org_info$pkg))

cat(sprintf("[ENRICH] Loaded %s\n", org_info$pkg))

# -----------------------------------------------------------------------------
# Read DA results
# -----------------------------------------------------------------------------
cat("[ENRICH] Reading DA results:", opt$da_results, "\n")
res <- fread(opt$da_results, sep = "\t", header = TRUE, data.table = TRUE)

required_cols <- c("protein_id", "gene_names", "significant", "comparison")
missing <- setdiff(required_cols, colnames(res))
if (length(missing) > 0) {
  stop(sprintf("Missing columns in DA results: %s", paste(missing, collapse = ", ")))
}

comparisons <- unique(res$comparison)
cat(sprintf("[ENRICH] Comparisons: %s\n", paste(comparisons, collapse = ", ")))

# -----------------------------------------------------------------------------
# Helper: map identifiers to Entrez IDs
# -----------------------------------------------------------------------------
map_to_entrez <- function(gene_names_vec, protein_ids_vec, org_info, orgdb) {

  clean_genes <- gsub(";.*", "", gene_names_vec)
  clean_genes <- trimws(clean_genes)
  clean_genes <- unique(clean_genes[clean_genes != "" & !is.na(clean_genes)])

  mapped_entrez <- character(0)

  # Try primary ID type (gene symbol / ORF name)
  if (length(clean_genes) > 0) {
    mapped <- tryCatch(
      suppressMessages(suppressWarnings(
        bitr(clean_genes,
             fromType = org_info$id_type,
             toType   = "ENTREZID",
             OrgDb    = orgdb)
      )),
      error = function(e) NULL
    )
    if (!is.null(mapped) && nrow(mapped) > 0) {
      mapped_entrez <- unique(mapped$ENTREZID)
      cat(sprintf("[ENRICH]   Mapped via %s: %d / %d genes\n",
                  org_info$id_type, length(mapped_entrez), length(clean_genes)))
    }
  }

  # Fallback: try UniProt IDs if symbol mapping failed or returned few hits
  if (length(mapped_entrez) < 3) {
    uniprot_ids <- gsub(";.*", "", protein_ids_vec)
    uniprot_ids <- trimws(uniprot_ids)
    uniprot_ids <- unique(uniprot_ids[uniprot_ids != "" & !is.na(uniprot_ids)])

    # Only try UNIPROT fallback if it's a valid keytype for this OrgDb
    available_keytypes <- tryCatch(
      keytypes(orgdb), error = function(e) character(0)
    )

    if ("UNIPROT" %in% available_keytypes && length(uniprot_ids) > 0) {
      mapped_u <- tryCatch(
        suppressMessages(suppressWarnings(
          bitr(uniprot_ids,
               fromType = "UNIPROT",
               toType   = "ENTREZID",
               OrgDb    = orgdb)
        )),
        error = function(e) NULL
      )
      if (!is.null(mapped_u) && nrow(mapped_u) > 0) {
        new_entrez <- unique(mapped_u$ENTREZID)
        cat(sprintf("[ENRICH]   Fallback via UNIPROT: %d IDs\n", length(new_entrez)))
        mapped_entrez <- unique(c(mapped_entrez, new_entrez))
      }
    }
  }

  mapped_entrez
}

# -----------------------------------------------------------------------------
# Helper: run GO + KEGG for one gene set
# -----------------------------------------------------------------------------
run_enrichment <- function(entrez_ids, universe_entrez, label,
                           output_dir, org_info, orgdb,
                           fdr_cutoff, min_genes) {

  if (length(entrez_ids) < min_genes) {
    cat(sprintf("[ENRICH]   Skipping %s: only %d genes mapped (need >=%d)\n",
                label, length(entrez_ids), min_genes))
    return(NULL)
  }

  results <- list()

  # GO Biological Process
  cat(sprintf("[ENRICH]   GO BP for %s (%d genes)\n", label, length(entrez_ids)))
  go_result <- tryCatch(
    suppressMessages(enrichGO(
      gene          = entrez_ids,
      universe      = universe_entrez,
      OrgDb         = orgdb,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = fdr_cutoff,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )),
    error = function(e) {
      cat(sprintf("[ENRICH]   GO BP failed: %s\n", e$message)); NULL
    }
  )

  if (!is.null(go_result) && nrow(go_result@result) > 0) {
    sig_go <- go_result@result[go_result@result$p.adjust < fdr_cutoff, ]
    if (nrow(sig_go) > 0) {
      go_df <- as.data.table(go_result@result)
      go_df[, comparison := label]
      go_df[, analysis   := "GO_BP"]
      results[["go"]] <- go_df

      n_show <- min(20, nrow(sig_go))
      p <- dotplot(go_result, showCategory = n_show, font.size = 9) +
        labs(title = sprintf("GO BP: %s", gsub("_", " ", label))) +
        theme(plot.title = element_text(size = 11))
      plot_path <- file.path(output_dir, sprintf("go_dotplot_%s.pdf", label))
      pdf(plot_path, width = 10, height = max(5, n_show * 0.4 + 2))
      print(p)
      dev.off()
      cat(sprintf("[ENRICH]   Wrote: %s\n", plot_path))
    } else {
      cat(sprintf("[ENRICH]   No significant GO terms for %s\n", label))
    }
  }

  # KEGG
  cat(sprintf("[ENRICH]   KEGG for %s\n", label))
  kegg_result <- tryCatch(
    suppressMessages(enrichKEGG(
      gene          = entrez_ids,
      universe      = universe_entrez,
      organism      = org_info$kegg,
      pAdjustMethod = "BH",
      pvalueCutoff  = fdr_cutoff,
      qvalueCutoff  = 0.2
    )),
    error = function(e) {
      cat(sprintf("[ENRICH]   KEGG failed: %s\n", e$message)); NULL
    }
  )

  if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
    sig_kegg <- kegg_result@result[kegg_result@result$p.adjust < fdr_cutoff, ]
    if (nrow(sig_kegg) > 0) {
      kegg_df <- as.data.table(kegg_result@result)
      kegg_df[, comparison := label]
      kegg_df[, analysis   := "KEGG"]
      results[["kegg"]] <- kegg_df

      n_show <- min(20, nrow(sig_kegg))
      p <- dotplot(kegg_result, showCategory = n_show, font.size = 9) +
        labs(title = sprintf("KEGG: %s", gsub("_", " ", label))) +
        theme(plot.title = element_text(size = 11))
      plot_path <- file.path(output_dir, sprintf("kegg_dotplot_%s.pdf", label))
      pdf(plot_path, width = 10, height = max(5, n_show * 0.4 + 2))
      print(p)
      dev.off()
      cat(sprintf("[ENRICH]   Wrote: %s\n", plot_path))
    } else {
      cat(sprintf("[ENRICH]   No significant KEGG pathways for %s\n", label))
    }
  }

  if (length(results) == 0) return(NULL)
  rbindlist(results, use.names = TRUE, fill = TRUE)
}

# -----------------------------------------------------------------------------
# Build universe
# -----------------------------------------------------------------------------
cat("[ENRICH] Building background gene universe\n")
universe_entrez <- map_to_entrez(res$gene_names, res$protein_id,
                                  org_info, orgdb)
cat(sprintf("[ENRICH] Universe: %d Entrez IDs\n", length(universe_entrez)))

out_path <- file.path(opt$output_dir, "enrichment_results.tsv")

if (length(universe_entrez) == 0) {
  cat("[ENRICH] WARNING: No genes mapped to Entrez IDs.\n")
  cat(sprintf("[ENRICH] Check that --organism %s matches your data.\n", opt$organism))
  empty <- data.table(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    comparison = character(), analysis = character()
  )
  fwrite(empty, out_path, sep = "\t", quote = FALSE)
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Run per comparison
# -----------------------------------------------------------------------------
all_enrichment <- list()

for (comp in comparisons) {
  cat(sprintf("[ENRICH] Processing: %s\n", comp))
  sig <- res[comparison == comp & significant == TRUE]
  cat(sprintf("[ENRICH]   Significant proteins: %d\n", nrow(sig)))

  if (nrow(sig) == 0) next

  entrez_ids <- map_to_entrez(sig$gene_names, sig$protein_id,
                               org_info, orgdb)
  cat(sprintf("[ENRICH]   Mapped to Entrez: %d\n", length(entrez_ids)))

  enrich_res <- run_enrichment(
    entrez_ids      = entrez_ids,
    universe_entrez = universe_entrez,
    label           = comp,
    output_dir      = opt$output_dir,
    org_info        = org_info,
    orgdb           = orgdb,
    fdr_cutoff      = opt$fdr_cutoff,
    min_genes       = opt$min_genes
  )

  if (!is.null(enrich_res)) all_enrichment[[comp]] <- enrich_res
}

# -----------------------------------------------------------------------------
# Write combined results
# -----------------------------------------------------------------------------
if (length(all_enrichment) > 0) {
  combined <- rbindlist(all_enrichment, use.names = TRUE, fill = TRUE)
  fwrite(combined, out_path, sep = "\t", quote = FALSE, na = "NA")
  cat(sprintf("[ENRICH] Wrote %d enrichment rows: %s\n", nrow(combined), out_path))
} else {
  cat("[ENRICH] No significant enrichment terms found.\n")
  empty <- data.table(
    ID = character(), Description = character(), GeneRatio = character(),
    BgRatio = character(), pvalue = numeric(), p.adjust = numeric(),
    qvalue = numeric(), geneID = character(), Count = integer(),
    comparison = character(), analysis = character()
  )
  fwrite(empty, out_path, sep = "\t", quote = FALSE)
}

cat("[ENRICH] Done.\n")
