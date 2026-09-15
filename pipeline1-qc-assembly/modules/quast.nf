process QUAST {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}/07_quast/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/quast-5.3.0.sif"

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
