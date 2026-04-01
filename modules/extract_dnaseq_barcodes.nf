#!/usr/bin/env nextflow

// =============================================================================
// DNA-seq clone barcode extraction module
// Supports two modes:
// 1. Whitelist mode (default): Use known barcode reference
// 2. Discovery mode: Two-pass approach to discover barcodes from data
// =============================================================================

process dnaseq_trim_reads {
    label 'medium'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
        path fastq_file

    output:
        path "${out_file}"

    script:
    out_file = fastq_file.getSimpleName() + ".${params.barcode_length}bp_5prime.fq.gz"
    
    """
    trim_galore --hardtrim5 ${params.barcode_length} \
                --cores ${task.cpus} $fastq_file \
                --gzip
    
    """
}

process dnaseq_filter_reads {
    label 'medium'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
        path fastq_file

    output:
        path "${out_file}"

    script:
    out_file = fastq_file.getSimpleName() + "_filtered.fq.gz"
    
    """
    fastp -i $fastq_file \
            -o $out_file \
            -q ${params.fastp_phred_for_qualified_reads} \
            -u ${params.fastp_percent_bases_unqualified} \
            -w ${task.cpus}
    """
}

process dnaseq_count_reads {
    // Add dummy adapters, run flexiplex discovery
    label 'small'

    input:
        path fastq_file

    output:
        path "${sample_name}_barcodes_counts.txt"

    script:
    sample_name = fastq_file.getSimpleName()
    fastq_w_adapter = sample_name + "_wDummyAdaptor.fq"
    """
    zcat $fastq_file | sed 's/^/START/g' | sed 's/START@/@/g' > ${fastq_w_adapter}
    
    flexiplex \
        -x "START" \
        -b ${params.barcode_length_chr} \
        -u "" \
        -x "" \
        -f 0 \
        -n $sample_name \
        -p ${task.cpus} \
        ${fastq_w_adapter}
    
    """
}

process dnaseq_split_reads_to_chunks {
    // break up the barcodes into chunks
    label 'small'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
        path barcode_counts

    output:
        path "${outdir}/${barcode_counts.baseName}_chunk*.fasta"

    script:
        outdir = "${barcode_counts.baseName}_unmapped_chunks"

    """
    mkdir ${outdir}
    dnaseq_split_reads.py --barcode_file ${barcode_counts} \
                                --sample_name ${barcode_counts.baseName} \
                                --n_chunks ${params.n_chunks} \
                                --outdir ${outdir}
    
    """
}

// =============================================================================
// Discovery mode processes for DNA-seq
// =============================================================================

process dnaseq_discover_barcodes {
    // Pass 1: Run flexiplex WITHOUT -k to discover barcodes from data
    // Uses -f 0 for strict flanking sequence match
    label 'small'

    input:
        path fastq_file

    output:
        path "${sample_name}_barcodes_counts.txt"

    script:
    sample_name = fastq_file.getSimpleName()
    fastq_w_adapter = sample_name + "_wDummyAdaptor.fq"
    """
    zcat $fastq_file | sed 's/^/START/g' | sed 's/START@/@/g' > ${fastq_w_adapter}
    
    # Run flexiplex in discovery mode (no -k flag)
    flexiplex \
        -x "START" \
        -b ${params.barcode_length_chr} \
        -u "" \
        -x "" \
        -f 0 \
        -n $sample_name \
        -p ${task.cpus} \
        ${fastq_w_adapter}
    
    """
}

process dnaseq_filter_discovered_barcodes {
    // Optionally filter discovered barcodes using flexiplex-filter knee-plot method
    // When params.filter_discovered_barcodes = false (default), all discovered
    // barcodes are kept using --no-inflection.
    // When params.filter_discovered_barcodes = true, knee-plot filtering is applied.
    label 'small'
    
    input:
        path barcode_counts

    output:
        path "filtered_barcodes.txt"

    """
    #!/usr/bin/bash
    
    # Run flexiplex-filter:
    # - filter_discovered_barcodes = false: --no-inflection keeps ALL discovered barcodes
    # - filter_discovered_barcodes = true:  knee-plot filtering removes low-count barcodes
    flexiplex-filter \
        ${params.filter_discovered_barcodes ? '' : '--no-inflection'} \
        --outfile filtered_barcodes.txt \
        ${barcode_counts}
    """
}

// =============================================================================
// Mapping processes (whitelist and discovery mode)
// =============================================================================

process dnaseq_map_barcodes {
    // Map barcodes using known reference (whitelist mode)
    // Ran flexiplex per fasta chunk
    // Then combine the counting of read (flexiplex discovery)
    // and the mapped barcode
    label "${params.mapping_process_profile}"
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
        path unmapped_fasta

    output:
        path "${out_file}"

    script:
    sample_name = unmapped_fasta.getSimpleName()
    mapped_chunk = sample_name + "_reads_barcodes.txt"
    out_file = sample_name + "_mapped.csv"
    """

    flexiplex \
        -x "START" \
        -b ${params.barcode_length_chr} \
        -u "" \
        -x "" \
        -f 0 \
        -n ${sample_name} \
        -k ${params.clone_barcodes_reference} \
        -e ${params.barcode_edit_distance} \
        -p ${task.cpus} \
        ${unmapped_fasta}

    dnaseq_combine_read_cnt_map.py --unmapped_chunk ${unmapped_fasta} \
                                --mapped_chunk ${mapped_chunk} \
                                --out_file ${out_file}
    """
}

process dnaseq_map_with_discovered_barcodes {
    // Pass 2: Map barcodes using discovered/filtered barcode list
    label "${params.mapping_process_profile}"
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
        path unmapped_fasta
        path discovered_barcodes

    output:
        path "${out_file}"

    script:
    sample_name = unmapped_fasta.getSimpleName()
    mapped_chunk = sample_name + "_reads_barcodes.txt"
    out_file = sample_name + "_mapped.csv"
    """

    flexiplex \
        -x "START" \
        -b ${params.barcode_length_chr} \
        -u "" \
        -x "" \
        -f 0 \
        -n ${sample_name} \
        -k ${discovered_barcodes} \
        -e ${params.barcode_edit_distance} \
        -p ${task.cpus} \
        ${unmapped_fasta}

    dnaseq_combine_read_cnt_map.py --unmapped_chunk ${unmapped_fasta} \
                                --mapped_chunk ${mapped_chunk} \
                                --out_file ${out_file}
    """
}

process dnaseq_collapse_barcodes {
    label 'small'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
        path mapped_reads

    output:
        path "clone_barcode_counts.csv"

    script:
    """
    dnaseq_count_barcodes.py . ${mapped_reads}
    """
}
