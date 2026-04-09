#!/bin/bash nextflow

// =============================================================================
// NextClone - Clonal barcode extraction pipeline
// Supports both DNAseq and scRNAseq modes
// 
// Two barcode identification approaches:
// 1. Whitelist mode (default): Use known barcode reference (clone_barcodes_reference)
// 2. Discovery mode: Two-pass approach to discover barcodes from data
//    - Pass 1: Run Flexiplex without -k to discover barcodes (-f 0 for strict match)
//    - Filter: Use flexiplex-filter (knee-plot method)
//    - Pass 2: Run Flexiplex with the discovered/filtered barcode list
// =============================================================================

params.barcode_length_chr = '?' * params.barcode_length

// Import DNAseq processes
include { 
    dnaseq_trim_reads;
    dnaseq_filter_reads;
    dnaseq_count_reads;
    dnaseq_split_reads_to_chunks;
    dnaseq_map_barcodes;
    dnaseq_discover_barcodes;
    dnaseq_filter_discovered_barcodes;
    dnaseq_map_with_discovered_barcodes;
    dnaseq_collapse_barcodes
} from "./modules/extract_dnaseq_barcodes"

// Import scRNAseq processes
include { 
    sc_get_unmapped_reads;
    sc_remove_low_qual_reads;
    sc_retain_reads_with_CB_tag;
    sc_split_unmapped_reads;
    sc_map_unmapped_reads;
    sc_discover_barcodes;
    sc_merge_discovered_barcodes;
    sc_map_with_discovered_barcodes;
    sc_merge_barcodes;
    generate_report;
    generate_run_log
} from "./modules/extract_sc_clone_barcodes"

workflow {

    // =============================================================================
    // Parameter validation
    // =============================================================================
    
    // Validate: discovery_mode = false requires clone_barcodes_reference
    if (!params.discovery_mode && !params.clone_barcodes_reference) {
        error """
        ERROR: Parameter 'clone_barcodes_reference' is required when 'discovery_mode = false'.
        
        Either:
        1. Provide a barcode whitelist: --clone_barcodes_reference /path/to/barcodes.txt
        2. Enable discovery mode: --discovery_mode true
        
        See documentation for details: https://phipsonlab.github.io/NextClone/
        """
    }
    
    // Validate: discovery_mode = true should not use clone_barcodes_reference (warn if provided)
    if (params.discovery_mode && params.clone_barcodes_reference) {
        log.warn """
        WARNING: 'clone_barcodes_reference' is ignored when 'discovery_mode = true'.
        Barcodes will be discovered from the data instead.
        """
    }



    if (params.mode == 'DNAseq') {
        
        if (params.discovery_mode) {
            // =========================================
            // Discovery mode workflow for DNAseq
            // =========================================
            
            // Preprocessing: trim and filter reads
            ch_filtered_reads = Channel.fromPath("${params.dnaseq_fastq_files}/*.fastq.gz") | 
                dnaseq_trim_reads |
                dnaseq_filter_reads
            
            // Pass 1: Discover barcodes from filtered reads
            ch_discovered = dnaseq_discover_barcodes(ch_filtered_reads)
            
            // Combine all discovered barcode counts and filter using knee-plot method
            ch_filtered_barcodes = dnaseq_filter_discovered_barcodes(
                ch_discovered.collectFile(name: 'combined_barcodes_counts.txt')
            )
            
            // Pass 2: Re-read files, preprocess, split, and map with discovered barcodes
            ch_barcode_chunks = Channel.fromPath("${params.dnaseq_fastq_files}/*.fastq.gz") |
                dnaseq_count_reads |
                dnaseq_split_reads_to_chunks
            
            ch_barcode_mappings = dnaseq_map_with_discovered_barcodes(
                ch_barcode_chunks.flatten(),
                ch_filtered_barcodes.first()
            )
            
            dnaseq_collapse_barcodes(ch_barcode_mappings.collect())
            
        } else {
            // =========================================
            // Whitelist mode workflow (original behavior)
            // =========================================
            
            ch_barcode_chunks = Channel.fromPath("${params.dnaseq_fastq_files}/*.fastq.gz") | 
                dnaseq_trim_reads |
                dnaseq_filter_reads |
                dnaseq_count_reads |
                dnaseq_split_reads_to_chunks
            
            ch_barcode_mappings = dnaseq_map_barcodes(ch_barcode_chunks.flatten())
            dnaseq_collapse_barcodes(ch_barcode_mappings.collect())
        }

    } 
    
    if (params.mode == 'scRNAseq') {
        
        // Initial preprocessing: get unmapped reads with cell barcodes
        ch_unmapped_fastas = Channel.fromPath("${params.scrnaseq_bam_files}/*.bam") | 
            sc_get_unmapped_reads |
            sc_remove_low_qual_reads |
            sc_retain_reads_with_CB_tag |
            sc_split_unmapped_reads
        
        if (params.discovery_mode) {
            // =========================================
            // Discovery mode workflow for scRNAseq
            // =========================================
            
            // Pass 1: Discover barcodes from each chunk
            ch_discovered = sc_discover_barcodes(ch_unmapped_fastas[0].flatten())
            
            // Combine and optionally filter discovered barcodes
            // sc_merge_discovered_barcodes handles both cases via params.filter_discovered_barcodes:
            // - false (default): --no-inflection keeps ALL discovered barcodes
            // - true: knee-plot filtering removes low-count barcodes
            ch_filtered_barcodes = sc_merge_discovered_barcodes(
                ch_discovered.collect()
            )
            
            // Pass 2: Map reads using discovered barcode list
            ch_mapped_fastas = sc_map_with_discovered_barcodes(
                ch_unmapped_fastas[0].flatten(),
                ch_filtered_barcodes.first()
            )
            
            ch_clone_barcodes = sc_merge_barcodes(ch_mapped_fastas.collect())
            generate_report(ch_clone_barcodes)
            generate_run_log(ch_clone_barcodes)
            
        } else {
            // =========================================
            // Whitelist mode workflow (original behavior)
            // =========================================
            
            ch_mapped_fastas = sc_map_unmapped_reads(ch_unmapped_fastas[0].flatten())
            ch_clone_barcodes = sc_merge_barcodes(ch_mapped_fastas.collect())
            generate_report(ch_clone_barcodes)
            generate_run_log(ch_clone_barcodes)
        }
    }
}
