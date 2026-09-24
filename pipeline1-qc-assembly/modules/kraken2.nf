process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/05_kraken2/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/Kraken/2.1.2/Kraken-2.1.2.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path db

    output:
    tuple val(meta), path("${meta.id}.kraken2.report.txt"), emit: report
    path "${meta.id}.kraken2.classified.txt", emit: classification

    script:
    """
    kraken2 --db ${db} \\
        --threads ${task.cpus} \\
        --confidence ${params.kraken2_confidence} \\
        --paired ${reads1} ${reads2} \\
        --report ${meta.id}.kraken2.report.txt \\
        --output ${meta.id}.kraken2.classified.txt \\
        --gzip-compressed
    """
}

// Re-estimates abundance from the Kraken2 report. Needs
// database<read_length>mers.kmer_distrib inside the Kraken2 database.
process BRACKEN {
    tag "$meta.id"
    label 'process_low'
    publishDir { "${params.outdir}/05_kraken2/${meta.id}" }, mode: 'copy'

    container params.bracken_sif

    input:
    tuple val(meta), path(kraken2_report)
    path db

    output:
    tuple val(meta), path("${meta.id}.bracken_${params.bracken_level}.tsv"), emit: abundance
    path "${meta.id}.bracken_${params.bracken_level}.report.txt", emit: report

    script:
    """
    bracken \\
        -d ${db} \\
        -i ${kraken2_report} \\
        -r ${params.bracken_read_length} \\
        -l ${params.bracken_level} \\
        -o ${meta.id}.bracken_${params.bracken_level}.tsv \\
        -w ${meta.id}.bracken_${params.bracken_level}.report.txt
    """
}
