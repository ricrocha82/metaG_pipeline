process FASTQC_RAW {
    tag "$meta.id"
    label 'process_low'
    publishDir { "${params.outdir}/01_fastqc_raw/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/FastQC/0.11.8/FastQC-0.11.8.sif"

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
    publishDir { "${params.outdir}/04_fastqc_trimmed/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/FastQC/0.11.8/FastQC-0.11.8.sif"

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
