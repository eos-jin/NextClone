#!/usr/bin/env nextflow

// =============================================================================
// LARRY barcode extraction module
// LARRY (Lineage And RNA Recovery) is a Cas9-based lineage tracing method
// from the Klein Lab (Weinreb et al., Science 2020)
//
// This module extracts LARRY barcodes from paired R1/R2 FASTQ files:
// - R1: Cell barcode (16bp) + UMI (8bp)
// - R2: LARRY barcode (40bp after prefix GTTGCTAGGAGAGACCATATG)
//
// Output is compatible with NextClone's DNAseq workflow for downstream
// barcode collapsing, filtering, and clone assignment.
// =============================================================================

process larry_extract_barcodes {
    // Extract LARRY barcodes from paired R1/R2 FASTQ files
    label 'medium'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
    path r1_fastq
    path r2_fastq

    output:
    path "${sample_name}_larry_barcodes.fastq.gz"
    path "${sample_name}_larry_stats.txt"

    script:
    sample_name = r1_fastq.baseName.replaceAll(/_R1.*$/, '')
    
    """
    larry_extract_barcodes.py \\
        --r1 $r1_fastq \\
        --r2 $r2_fastq \\
        --output ${sample_name}_larry_barcodes.fastq.gz \\
        --prefix ${params.larry_prefix} \\
        --cell-bc-len ${params.larry_cell_bc_len} \\
        --umi-len ${params.larry_umi_len} \\
        --larry-bc-len ${params.larry_bc_len} \\
        2> ${sample_name}_larry_stats.txt
    """
}

process larry_filter_and_cluster {
    // Complete LARRY filtering and clustering workflow
    // Implements the full pipeline from Weinreb et al. Science 2020:
    // 1. Count (cell_bc, umi, larry_bc) tuples
    // 2. Filter by minimum reads per tuple
    // 3. Cluster LARRY barcodes by Hamming distance
    // 4. Count UMIs per (cell, barcode) combination
    // 5. Filter by minimum UMIs per (cell, barcode)
    // 6. Output clone assignments per cell
    label 'medium'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
    path larry_fastq
    val sample_id

    output:
    path "${sample_id}_larry_clones.csv"
    path "${sample_id}_larry_filter_stats.txt"

    script:
    """
    larry_filter_and_cluster.py \
        ${larry_fastq} \
        ${sample_id}_larry_clones.csv \
        --sample ${sample_id} \
        --min-reads ${params.larry_min_reads} \
        --min-umis ${params.larry_min_umis} \
        --max-hamming ${params.larry_max_hamming} \
        2> ${sample_id}_larry_filter_stats.txt
    """
}
