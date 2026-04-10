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
    // barcodes are kept (no filtering). This is recommended for lineage tracing 
    // where singleton clones are biologically meaningful and should not be discarded.
    // When params.filter_discovered_barcodes = true, the knee-plot inflection point
    // method is used to remove low-count/noisy barcodes.
    label 'small'
    
    input:
        path barcode_counts_files

    output:
        path "all_barcodes.txt"
        path "filtered_barcodes.txt"

    """
    #!/usr/bin/bash
    
    # Combine all barcode counts files
    # Sum counts for same barcodes across chunks
    cat ${barcode_counts_files} | \
        awk '{counts[\$1] += \$2} END {for (bc in counts) print bc "\\t" counts[bc]}' | \
        sort -k2 -nr > combined_barcodes_counts.txt
    
    # Save ALL discovered barcodes (no filtering) - useful for debugging and QC
    # Header: #barcode = lineage tracing barcode sequence, count = number of reads supporting this barcode
    echo -e "#barcode\\tcount" > all_barcodes.txt
    echo "# barcode: lineage tracing barcode sequence" >> all_barcodes.txt
    echo "# count: number of reads supporting this barcode" >> all_barcodes.txt
    cat combined_barcodes_counts.txt >> all_barcodes.txt
    
    # IMPORTANT: flexiplex-filter has default bounds (min-rank=50, max-rank=95th percentile)
    # Even with --no-inflection, it still filters! So we must NOT call it when filtering is disabled.
    # See: https://davidsongroup.github.io/flexiplex/tutorial.html
    
    if [ "${params.filter_discovered_barcodes}" = "true" ]; then
        # Run knee-plot inflection point filtering
        flexiplex-filter \
            --outfile filtered_barcodes.txt.tmp \
            combined_barcodes_counts.txt
        # Add header with explanation
        echo "#barcode\tcount" > filtered_barcodes.txt
        echo "# barcode: lineage tracing barcode sequence" >> filtered_barcodes.txt
        echo "# count: number of reads supporting this barcode" >> filtered_barcodes.txt
        cat filtered_barcodes.txt.tmp >> filtered_barcodes.txt
        rm -f filtered_barcodes.txt.tmp
    else
        # NO filtering at all - just copy the combined file directly
        # This preserves ALL barcodes including singletons
        echo "#barcode\tcount" > filtered_barcodes.txt
        echo "# barcode: lineage tracing barcode sequence" >> filtered_barcodes.txt
        echo "# count: number of reads supporting this barcode" >> filtered_barcodes.txt
        cat combined_barcodes_counts.txt >> filtered_barcodes.txt
    fi
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

process generate_run_log {
    // Generate run log with parameters and command line
    // Saved to publish_dir for reproducibility
    label 'small'
    
    publishDir params.publish_dir, mode: params.publish_dir_mode
    
    input:
        path clone_barcodes
    
    output:
        path "run_log.txt"
    
    script:
        timestamp = new Date().format('yyyy-MM-dd HH:mm:ss')
    """
    # Get software versions
    NF_VERSION=\$(nextflow -version 2>&1 | head -1 || echo "unknown")
    FLEXIPLEX_VERSION=\$(flexiplex --version 2>&1 | head -1 || echo "unknown")
    PYTHON_VERSION=\$(python3 --version 2>&1 || echo "unknown")
    
    # Get git info if available
    GIT_COMMIT=\$(git rev-parse HEAD 2>/dev/null || echo "Not a git repo")
    GIT_BRANCH=\$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
    
    cat > run_log.txt << EOF
# NextClone Run Log
# Generated: ${timestamp}

## Software Versions
Nextflow: \${NF_VERSION}
Flexiplex: \${FLEXIPLEX_VERSION}
Python: \${PYTHON_VERSION}

## Code Version
Git commit: \${GIT_COMMIT}
Git branch: \${GIT_BRANCH}

## Command
nextflow run ${projectDir}/main.nf \\
    --mode ${params.mode} \\
    --discovery_mode ${params.discovery_mode} \\
    --filter_discovered_barcodes ${params.filter_discovered_barcodes} \\
    --barcode_edit_distance ${params.barcode_edit_distance} \\
    --adapter_edit_distance ${params.adapter_edit_distance} \\
    --n_chunks ${params.n_chunks} \\
    --publish_dir ${params.publish_dir}

## Parameters
mode = ${params.mode}
discovery_mode = ${params.discovery_mode}
filter_discovered_barcodes = ${params.filter_discovered_barcodes}
barcode_edit_distance = ${params.barcode_edit_distance}
adapter_edit_distance = ${params.adapter_edit_distance}
barcode_length = ${params.barcode_length}
n_chunks = ${params.n_chunks}
publish_dir = ${params.publish_dir}

## Output Files
- all_barcodes.txt: All discovered barcodes (no filtering)
- filtered_barcodes.txt: Barcodes after filtering (same as all_barcodes.txt if filter_discovered_barcodes=false)
- clone_barcodes.csv: Final clone assignments to cells
- nextclone_qc_report.html: Interactive QC dashboard

## Notes
- all_barcodes.txt contains ALL barcodes discovered in Pass 1, including singletons
- filtered_barcodes.txt applies knee-plot filtering only if filter_discovered_barcodes=true
- For lineage tracing, we recommend filter_discovered_barcodes=false to retain rare clones
EOF
    """
}

process generate_report {
    // Generate interactive HTML dashboard from clone_barcodes.csv
    // Uses reports/generate_report.py (pure Python stdlib, no pip installs)
    label 'small'

    publishDir params.publish_dir, mode: params.publish_dir_mode

    input:
        path clone_barcodes

    output:
        path "nextclone_qc_report.html"

    script:
        title = params.report_title ?: "NextClone QC Report — ${new Date().format('yyyy-MM-dd')}"
    """
    python3 ${projectDir}/reports/generate_report.py \
        ${clone_barcodes} \
        --output nextclone_qc_report.html \
        --title "${title}"
    """
}
