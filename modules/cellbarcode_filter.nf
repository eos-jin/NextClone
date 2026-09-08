#!/usr/bin/env nextflow

// =============================================================================
// CellBarcode filtering module for NextClone
// Uses the official CellBarcode R package (Sun et al. 2024)
// https://github.com/wenjie1991/CellBarcode
// =============================================================================

process cellbarcode_filter {
    label 'small'
    conda "${projectDir}/conda_env/cellbarcode_env.yaml"

    input:
    path barcode_counts

    output:
    path "filtered_barcodes.txt"
    path "cellbarcode_stats.txt"

    script:
    def method = params.cellbarcode_method ?: 'auto'
    def threshold = params.cellbarcode_threshold ?: 'NA'
    def cluster_dist = params.cellbarcode_cluster_distance ?: 1
    def min_count = params.cellbarcode_min_count ?: 1

    """
    cellbarcode_filter.R \\
        ${barcode_counts} \\
        filtered_barcodes.txt \\
        ${method} \\
        ${threshold} \\
        ${cluster_dist} \\
        ${min_count} \\
        2> cellbarcode_stats.txt
    """
}
