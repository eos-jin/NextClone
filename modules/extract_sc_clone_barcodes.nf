#!/usr/bin/env nextflow

// =============================================================================
// Single-cell RNA-seq clone barcode extraction module
// Supports two modes:
// 1. Whitelist mode (default): Use known barcode reference
// 2. Discovery mode: Two-pass approach to discover barcodes from data
// =============================================================================

process sc_get_unmapped_reads {
    // Using sambamba
    module 'sambamba'
    label 'medium_mem'

    input:
        path bam_file

    output:
        path "${out_bam_file}"

    script:
        out_bam_file = "${bam_file.baseName}_reads_unmapped.bam"
    """
    sambamba view -t ${task.cpus} -F "unmapped" -f bam -o ${out_bam_file} ${bam_file}
    """
}

process sc_remove_low_qual_reads {
    label 'small_mem'
    conda "${projectDir}/conda_env/extract_sc_env.yaml"
    
    input:
        path unmapped_bam

    output:
        path "${out_bam_file}"

    script:
        out_bam_file = "${unmapped_bam.baseName}_filtered.bam"

    """
    sc_remove_low_qual_reads.py ${unmapped_bam} ${params.phred_thres} ${out_bam_file}
    """
}

process sc_retain_reads_with_CB_tag {
    // Using sambamba
    module 'sambamba'
    label 'medium_mem'

    input:
        path bam_file

    output:
        path "${out_bam_file}"

    script:
        out_bam_file = "${bam_file.baseName}_withCB.bam"
    
    """
    #!/usr/bin/bash
    sambamba view \
        -F "([CB] != null and [UB] != null)" \
        -t ${task.cpus} \
        -f bam \
        -o ${out_bam_file} \
        ${bam_file}
    """
}

process sc_split_unmapped_reads {
    label 'small_mem'
    conda "${projectDir}/conda_env/extract_sc_env.yaml"

    input:
        path unmapped_bam 

    output:
        path "${outdir}/${unmapped_bam.baseName}_unmapped_chunk_*.fasta"

    script:
        outdir = "${unmapped_bam.baseName}_unmapped_chunks"
    """
    mkdir ${outdir}
    
    sc_split_reads.py \
        --input_bam_filename ${unmapped_bam} \
        --outdir ${outdir} \
        --n_chunks ${params.n_chunks}
    """
}

// =============================================================================
// Discovery mode processes (Pass 1 and filtering)
// =============================================================================

process sc_discover_barcodes {
    // Pass 1: Run flexiplex WITHOUT -k to discover barcodes from data
    // Uses -f 0 for strict flanking sequence match to reduce errors
    label "${params.mapping_process_profile}"
    
    input:
        path unmapped_fasta

    output:
        path "${unmapped_fasta.baseName}_barcodes_counts.txt"

    script:
    """
    #!/usr/bin/bash
    
    # Run flexiplex in discovery mode (no -k flag)
    # -f 0: strict flanking sequence match (reduces barcode errors)
    flexiplex \
        -x "${params.adapter_5prime}" \
        -b ${params.barcode_length_chr} \
        -u "" \
        -x "${params.adapter_3prime}" \
        -f 0 \
        -n ${unmapped_fasta.baseName} \
        -p ${task.cpus} \
        ${unmapped_fasta}
    """
}

process sc_filter_discovered_barcodes {
    // Filter discovered barcodes using flexiplex-filter
    // Uses knee-plot inflection point method
    // Optionally intersects with 10x whitelist if provided
    label 'small'
    
    input:
        path barcode_counts

    output:
        path "filtered_barcodes.txt"

    """
    #!/usr/bin/bash
    
    # Run flexiplex-filter to select high-quality barcodes
    # Uses knee-plot inflection point method
    flexiplex-filter \
        --outfile filtered_barcodes.txt \
        ${barcode_counts}
    """
}

process sc_merge_discovered_barcodes {
    // Merge barcode counts from all chunks and optionally filter using knee-plot
    // When params.filter_discovered_barcodes = false (default), all discovered
    // barcodes are kept using flexiplex-filter --no-inflection.
    // This is recommended for lineage tracing where singleton clones are biologically
    // meaningful and should not be discarded.
    // When params.filter_discovered_barcodes = true, the knee-plot inflection point
    // method is used to remove low-count/noisy barcodes.
    label 'small'
    
    input:
        path barcode_counts_files

    output:
        path "filtered_barcodes.txt"

    """
    #!/usr/bin/bash
    
    # Combine all barcode counts files
    # Sum counts for same barcodes across chunks
    cat ${barcode_counts_files} | \
        awk '{counts[\$1] += \$2} END {for (bc in counts) print bc "\\t" counts[bc]}' | \
        sort -k2 -nr > combined_barcodes_counts.txt
    
    # Run flexiplex-filter:
    # - filter_discovered_barcodes = false: --no-inflection keeps ALL discovered barcodes
    # - filter_discovered_barcodes = true:  knee-plot filtering removes low-count barcodes
    flexiplex-filter \
        ${params.filter_discovered_barcodes ? '' : '--no-inflection'} \
        --outfile filtered_barcodes.txt \
        combined_barcodes_counts.txt
    """
}

// =============================================================================
// Mapping processes (Pass 2 for discovery mode, or single pass for whitelist mode)
// =============================================================================

process sc_map_unmapped_reads {
    // Map reads to known barcode reference (whitelist mode)
    label "${params.mapping_process_profile}"

    input:
        path unmapped_fasta

    output:
        path "${unmapped_fasta.baseName}_reads_barcodes.txt"

    """
    #!/usr/bin/bash

    flexiplex \
            -x "${params.adapter_5prime}" \
            -b ${params.barcode_length_chr} \
            -u "" \
            -x "${params.adapter_3prime}" \
            -f ${params.adapter_edit_distance} \
            -e ${params.barcode_edit_distance} \
            -n ${unmapped_fasta.baseName} \
            -k ${params.clone_barcodes_reference} \
            -p ${task.cpus} \
            $unmapped_fasta
    
    """
}

process sc_map_with_discovered_barcodes {
    // Pass 2: Map reads using discovered/filtered barcode list
    label "${params.mapping_process_profile}"

    input:
        path unmapped_fasta
        path discovered_barcodes

    output:
        path "${unmapped_fasta.baseName}_reads_barcodes.txt"

    """
    #!/usr/bin/bash

    flexiplex \
            -x "${params.adapter_5prime}" \
            -b ${params.barcode_length_chr} \
            -u "" \
            -x "${params.adapter_3prime}" \
            -f ${params.adapter_edit_distance} \
            -e ${params.barcode_edit_distance} \
            -n ${unmapped_fasta.baseName} \
            -k ${discovered_barcodes} \
            -p ${task.cpus} \
            ${unmapped_fasta}
    
    """
}

process sc_merge_barcodes {
    label 'small_mem'
    conda "${projectDir}/conda_env/extract_sc_env.yaml"

    input:
        path mapped_reads

    output:
        path "${outfile}"

    script:
        outfile = "clone_barcodes.csv"

    """
    sc_merge_clone_barcodes.py ${mapped_reads} ${outfile}
    """
}
