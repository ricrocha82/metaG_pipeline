process FASTQC_RAW {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.qc_outdir}/01_fastqc_raw/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/FastQC-0.11.8.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads1} ${reads2}
    """
}

process FASTQC_TRIMMED {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.qc_outdir}/04_fastqc_trimmed/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/FastQC-0.11.8.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads1} ${reads2}
    """
}
