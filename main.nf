#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// RAPID — Reproducible Analysis Pipeline for Intensity-based
//         proteomics Data
// ============================================================
// Modules:
//   1. DOWNLOAD_DATA   — fetch MaxQuant proteinGroups.txt
//   2. QC_FILTER       — remove contaminants, reverse hits,
//                        low-coverage proteins
//   3. NORMALIZE       — log2 transform + median centering
//   4. PCA             — PCA plot of normalized matrix
//   5. DA_MULTI        — all pairwise DEqMS comparisons
//   6. ENRICHMENT      — GO/KEGG pathway enrichment
//   7. REPORT          — self-contained HTML report
//   8. VALIDATE        — output validation and integrity checks
// ============================================================


// ============================================================
// PROCESS 1: Download data
// ============================================================
process DOWNLOAD_DATA {
    tag "download"
    label 'process_low'

    publishDir "${params.outdir}/raw", mode: 'copy'

    output:
    path "proteinGroups.txt", emit: protein_groups

    script:
    """
    wget -q --tries=3 --timeout=60 \\
        "${params.data_url}" \\
        -O proteinGroups.txt

    [ -s proteinGroups.txt ] || { echo "ERROR: Download failed or empty"; exit 1; }

    head -1 proteinGroups.txt | grep -q "Protein IDs" || \\
        { echo "ERROR: Not a MaxQuant proteinGroups file"; exit 1; }

    echo "Download complete: \$(wc -l < proteinGroups.txt) lines"
    """
}

// ============================================================
// PROCESS 2: QC and filtering
// ============================================================
process QC_FILTER {
    tag "qc_filter"
    label 'process_low'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path protein_groups
    path sample_sheet

    output:
    path "filtered_matrix.tsv", emit: filtered_matrix
    path "qc_summary.txt",      emit: qc_summary

    script:
    """
    Rscript ${projectDir}/bin/01_qc_filter.R \\
        --input         ${protein_groups} \\
        --sample-sheet  ${sample_sheet} \\
        --min-valid     ${params.min_valid_values} \\
        --output        filtered_matrix.tsv \\
        --qc-summary    qc_summary.txt
    """
}

// ============================================================
// PROCESS 3: Normalization
// ============================================================
process NORMALIZE {
    tag "normalize"
    label 'process_low'

    publishDir "${params.outdir}/normalized", mode: 'copy'

    input:
    path filtered_matrix
    path sample_sheet

    output:
    path "normalized_matrix.tsv",     emit: normalized_matrix
    path "normalization_boxplot.png",  emit: boxplot

    script:
    """
    Rscript ${projectDir}/bin/02_normalize.R \\
        --input         ${filtered_matrix} \\
        --sample-sheet  ${sample_sheet} \\
        --output        normalized_matrix.tsv \\
        --boxplot       normalization_boxplot.png
    """
}

// ============================================================
// PROCESS 4: PCA
// ============================================================
process PCA {
    tag "pca"
    label 'process_low'

    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path normalized_matrix
    path sample_sheet

    output:
    path "pca_plot.png",    emit: pca_plot
    path "pca_coords.tsv",  emit: pca_coords

    script:
    """
    Rscript ${projectDir}/bin/03_pca.R \\
        --input         ${normalized_matrix} \\
        --sample-sheet  ${sample_sheet} \\
        --output-plot   pca_plot.png \\
        --output-data   pca_coords.tsv
    """
}

// ============================================================
// PROCESS 5: All pairwise differential abundance
// ============================================================
process DA_MULTI {
    tag "da_multi"
    label 'process_low'

    publishDir "${params.outdir}/results", mode: 'copy'

    input:
    path normalized_matrix
    path sample_sheet

    output:
    path "da_results_all_comparisons.tsv", emit: da_results
    path "volcano_*.png",                  emit: volcanos
    path "da_results_*_vs_*.tsv",          emit: da_per_comparison, optional: true

    script:
    """
    Rscript ${projectDir}/bin/04_da_multi.R \\
        --input         ${normalized_matrix} \\
        --sample-sheet  ${sample_sheet} \\
        --fdr-cutoff    ${params.fdr_cutoff} \\
        --fc-cutoff     ${params.fc_cutoff} \\
        --output-dir    .
    """
}

// ============================================================
// PROCESS 6: Pathway enrichment
// ============================================================
process ENRICHMENT {
    tag "enrichment"
    label 'process_low'

    publishDir "${params.outdir}/enrichment", mode: 'copy'

    input:
    path da_results

    output:
    path "enrichment_results.tsv",  emit: enrichment_results
    path "*.png",                   emit: enrichment_plots, optional: true

    script:
    """
    Rscript ${projectDir}/bin/05_enrichment.R \\
        --da-results    ${da_results} \\
        --output-dir    . \\
        --organism      ${params.organism} \\
        --fdr-cutoff    ${params.fdr_cutoff}
    """
}

// ============================================================
// PROCESS 7: HTML report
// ============================================================
process REPORT {
    tag "report"
    label 'process_low'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path qc_summary
    path boxplot
    path pca_plot
    path da_results
    path enrichment_results
    path sample_sheet

    output:
    path "proteomics_report.html", emit: report

    script:
    """
    Rscript ${projectDir}/bin/06_report.R \\
        --qc-summary           ${qc_summary} \\
        --boxplot              ${boxplot} \\
        --pca-plot             ${pca_plot} \\
        --da-results           ${da_results} \\
        --enrichment-results   ${enrichment_results} \\
        --sample-sheet         ${sample_sheet} \\
        --condition-a          "${params.condition_a}" \\
        --condition-b          "${params.condition_b}" \\
        --fdr-cutoff           ${params.fdr_cutoff} \\
        --fc-cutoff            ${params.fc_cutoff} \\
        --output               proteomics_report.html
    """
}

// ============================================================
// PROCESS 8: Output validation
// ============================================================
process VALIDATE {
    tag "validate"
    label 'process_low'

    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path filtered_matrix
    path normalized_matrix
    path pca_coords
    path da_results
    path enrichment_results
    path report
    path qc_summary
    path sample_sheet

    output:
    path "validation_report.txt", emit: validation_report

    script:
    """
    Rscript ${projectDir}/bin/07_validate.R \\
        --filtered-matrix    ${filtered_matrix} \\
        --normalized-matrix  ${normalized_matrix} \\
        --pca-coords         ${pca_coords} \\
        --da-results         ${da_results} \\
        --enrichment-results ${enrichment_results} \\
        --report             ${report} \\
        --qc-summary         ${qc_summary} \\
        --sample-sheet       ${sample_sheet} \\
        --output             validation_report.txt \\
        --min-proteins       ${params.min_valid_values}
    """
}

// ============================================================
// WORKFLOW
// ============================================================
workflow {

log.info """
=========================================
 R A P I D   P I P E L I N E
=========================================
Dataset URL  : ${params.data_url}
Sample sheet : ${params.sample_sheet}
Min valid    : ${params.min_valid_values}
FDR cutoff   : ${params.fdr_cutoff}
FC cutoff    : ${params.fc_cutoff}
Output dir   : ${params.outdir}
=========================================
""".stripIndent()

    if (!params.data_url)    error "ERROR: --data_url is required"
    if (!params.sample_sheet) error "ERROR: --sample_sheet is required"

    sample_sheet_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)

    DOWNLOAD_DATA()

    QC_FILTER(
        DOWNLOAD_DATA.out.protein_groups,
        sample_sheet_ch
    )

    NORMALIZE(
        QC_FILTER.out.filtered_matrix,
        sample_sheet_ch
    )

    PCA(
        NORMALIZE.out.normalized_matrix,
        sample_sheet_ch
    )

    DA_MULTI(
        NORMALIZE.out.normalized_matrix,
        sample_sheet_ch
    )

    ENRICHMENT(
        DA_MULTI.out.da_results
    )

    REPORT(
        QC_FILTER.out.qc_summary,
        NORMALIZE.out.boxplot,
        PCA.out.pca_plot,
        DA_MULTI.out.da_results,
        ENRICHMENT.out.enrichment_results,
        sample_sheet_ch
    )

    VALIDATE(
        QC_FILTER.out.filtered_matrix,
        NORMALIZE.out.normalized_matrix,
        PCA.out.pca_coords,
        DA_MULTI.out.da_results,
        ENRICHMENT.out.enrichment_results,
        REPORT.out.report,
        QC_FILTER.out.qc_summary,
        sample_sheet_ch
    )
}
