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

process larry_count_barcodes {
    // Count unique LARRY barcodes and their frequencies
    // Output format: barcode\tcount (compatible with NextClone DNAseq workflow)
    label 'small'
    conda "${projectDir}/conda_env/extract_dnaseq_env.yaml"

    input:
    path larry_fastq

    output:
    path "${sample_name}_barcodes_counts.txt"

    script:
    sample_name = larry_fastq.baseName.replaceAll(/_larry_barcodes.*$/, '')
    
    """
    #!/usr/bin/env python3
    import gzip
    from collections import Counter
    
    # Count unique barcodes
    barcodes = Counter()
    
    with gzip.open('$larry_fastq', 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            if seq:
                barcodes[seq] += 1
    
    # Write counts (sorted by frequency, descending)
    with open('${sample_name}_barcodes_counts.txt', 'w') as out:
        for bc, count in sorted(barcodes.items(), key=lambda x: -x[1]):
            out.write(f'{bc}\\t{count}\\n')
    
    print(f'Counted {len(barcodes)} unique barcodes from {sum(barcodes.values())} total reads')
    """
}

process larry_split_reads_to_chunks {
    // Split LARRY barcodes into chunks for parallel mapping
    // Reuses the existing dnaseq_split_reads.py script
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
    dnaseq_split_reads.py --barcode_file ${barcode_counts} \\
                                --sample_name ${barcode_counts.baseName} \\
                                --n_chunks ${params.n_chunks} \\
                                --outdir ${outdir}
    """
}
