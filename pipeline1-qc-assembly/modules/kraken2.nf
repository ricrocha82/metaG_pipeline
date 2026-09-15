process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/05_kraken2/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/Kraken-2.17.1.sif"

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
        --paired ${reads1} ${reads2} \\
        --report ${meta.id}.kraken2.report.txt \\
        --output ${meta.id}.kraken2.classified.txt \\
        --gzip-compressed
    """
}
