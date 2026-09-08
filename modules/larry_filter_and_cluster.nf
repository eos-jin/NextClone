process LARRY_FILTER_AND_CLUSTER {
    tag "$sample_id"
    label 'process_medium'
    
    conda "conda-forge::python=3.9"
    
    input:
    path fastq
    val sample_id
    
    output:
    path "${sample_id}_larry_clones.csv", emit: clones
    path "versions.yml", emit: versions
    
    script:
    """
    larry_filter_and_cluster.py \\
        $fastq \\
        ${sample_id}_larry_clones.csv \\
        --sample $sample_id \\
        --min-reads $params.larry_min_reads \\
        --min-umis $params.larry_min_umis \\
        --max-hamming $params.larry_max_hamming
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
