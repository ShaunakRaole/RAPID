#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ============================================================
// ProteomicsDA — Reproducible Proteomics Differential Abundance
// ============================================================
// Modules:
//   1. DOWNLOAD_DATA   — fetch MaxQuant proteinGroups.txt from PRIDE
//   2. QC_FILTER       — remove contaminants, reverse hits, high-missingness proteins
//   3. NORMALIZE       — log2 transform + median centering
//   4. DA_ANALYSIS     — DEqMS differential abundance analysis
//   5. REPORT          — HTML report with volcano plot and summary table
// ============================================================

log.info """
=========================================
 P R O T E O M I C S - D A   P I P E L I N E
=========================================
Dataset URL  : ${params.data_url}
Sample sheet : ${params.sample_sheet}
Condition A  : ${params.condition_a}
Condition B  : ${params.condition_b}
Min valid    : ${params.min_valid_values}
FDR cutoff   : ${params.fdr_cutoff}
FC cutoff    : ${params.fc_cutoff}
Output dir   : ${params.outdir}
=========================================
""".stripIndent()

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

    # Sanity check: file must exist and be non-empty
    [ -s proteinGroups.txt ] || { echo "ERROR: Download failed or file is empty"; exit 1; }

    # Sanity check: must look like a MaxQuant proteinGroups file
    head -1 proteinGroups.txt | grep -q "Protein IDs" || \\
        { echo "ERROR: File does not appear to be a MaxQuant proteinGroups.txt"; exit 1; }

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
    path "filtered_matrix.tsv",   emit: filtered_matrix
    path "qc_summary.txt",        emit: qc_summary

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
    path "normalized_matrix.tsv",      emit: normalized_matrix
    path "normalization_boxplot.pdf",  emit: boxplot

    script:
    """
    Rscript ${projectDir}/bin/02_normalize.R \\
        --input         ${filtered_matrix} \\
        --sample-sheet  ${sample_sheet} \\
        --output        normalized_matrix.tsv \\
        --boxplot       normalization_boxplot.pdf
    """
}

// ============================================================
// PROCESS 4: Differential abundance analysis
// ============================================================
process DA_ANALYSIS {
    tag "da_analysis"
    label 'process_low'

    publishDir "${params.outdir}/results", mode: 'copy'

    input:
    path normalized_matrix
    path sample_sheet

    output:
    path "da_results.tsv",     emit: da_results
    path "volcano_plot.pdf",   emit: volcano

    script:
    """
    Rscript ${projectDir}/bin/03_da_analysis.R \\
        --input         ${normalized_matrix} \\
        --sample-sheet  ${sample_sheet} \\
        --condition-a   "${params.condition_a}" \\
        --condition-b   "${params.condition_b}" \\
        --fdr-cutoff    ${params.fdr_cutoff} \\
        --fc-cutoff     ${params.fc_cutoff} \\
        --output        da_results.tsv \\
        --volcano       volcano_plot.pdf
    """
}

// ============================================================
// PROCESS 5: HTML report
// ============================================================
process REPORT {
    tag "report"
    label 'process_low'

    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path qc_summary
    path boxplot
    path da_results
    path volcano
    path sample_sheet

    output:
    path "proteomics_report.html", emit: report

    script:
    """
    Rscript ${projectDir}/bin/04_report.R \\
        --qc-summary    ${qc_summary} \\
        --boxplot       ${boxplot} \\
        --da-results    ${da_results} \\
        --volcano       ${volcano} \\
        --sample-sheet  ${sample_sheet} \\
        --condition-a   "${params.condition_a}" \\
        --condition-b   "${params.condition_b}" \\
        --fdr-cutoff    ${params.fdr_cutoff} \\
        --fc-cutoff     ${params.fc_cutoff} \\
        --output        proteomics_report.html
    """
}

// ============================================================
// WORKFLOW
// ============================================================
workflow {

    // Validate required parameters
    if (!params.data_url)      error "ERROR: --data_url is required"
    if (!params.sample_sheet)  error "ERROR: --sample_sheet is required"
    if (!params.condition_a)   error "ERROR: --condition_a is required"
    if (!params.condition_b)   error "ERROR: --condition_b is required"

    // Input channels
    sample_sheet_ch = Channel.fromPath(params.sample_sheet, checkIfExists: true)

    // Run pipeline
    DOWNLOAD_DATA()

    QC_FILTER(
        DOWNLOAD_DATA.out.protein_groups,
        sample_sheet_ch
    )

    NORMALIZE(
        QC_FILTER.out.filtered_matrix,
        sample_sheet_ch
    )

    DA_ANALYSIS(
        NORMALIZE.out.normalized_matrix,
        sample_sheet_ch
    )

    REPORT(
        QC_FILTER.out.qc_summary,
        NORMALIZE.out.boxplot,
        DA_ANALYSIS.out.da_results,
        DA_ANALYSIS.out.volcano,
        sample_sheet_ch
    )
}
