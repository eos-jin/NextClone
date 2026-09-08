#!/usr/bin/env nextflow

// =============================================================================
// CellBarcode filtering module for NextClone
// Implements filtering strategies from CellBarcode paper (Sun et al. 2024)
// =============================================================================

process cellbarcode_filter {
    label 'small'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
    path barcode_counts

    output:
    path "filtered_barcodes.txt"
    path "cellbarcode_stats.txt"

    script:
    def method = params.cellbarcode_method ?: 'auto'
    def threshold = params.cellbarcode_threshold ?: ''
    def cluster_dist = params.cellbarcode_cluster_distance ?: 1
    def min_count = params.cellbarcode_min_count ?: 1

    """
    cellbarcode_filter.py \\
        ${barcode_counts} \\
        filtered_barcodes.txt \\
        --method ${method} \\
        ${threshold ? "--threshold " + threshold : ''} \\
        --cluster-distance ${cluster_dist} \\
        --min-count ${min_count} \\
        2> cellbarcode_stats.txt
    """
}
