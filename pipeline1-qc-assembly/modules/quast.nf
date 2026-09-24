process QUAST {
    tag "$meta.id"
    label 'process_medium'
    publishDir { "${params.outdir}/07_quast/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/quast/5.0.2/quast_5.0.2.sif"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_quast"), emit: results

    script:
    """
    quast.py ${assembly} \\
        -o ${meta.id}_quast \\
        -t ${task.cpus} \\
        --min-contig ${params.min_contig_size ?: 1000}
    """
}
