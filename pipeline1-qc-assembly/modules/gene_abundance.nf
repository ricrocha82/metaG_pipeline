process GENE_ABUNDANCE {
    // NOTE: no subread/featureCounts .sif found in the OSC PAS1117 module
    // library -- using conda since there's no local container to point at.
    // If your group later adds one, swap this `conda` line for a
    // `container "/fs/project/PAS1117/modules/singularity/<name>.sif"` line.
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.qc_outdir}/10_gene_abundance/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/subread-2.0.6.sif"

    input:
    tuple val(meta), path(gff), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.gene_abundance.tsv"), emit: abundance

    script:
    """
    featureCounts -T ${task.cpus} -p --countReadPairs \\
        -t CDS -g ID \\
        -a ${gff} -F GFF \\
        -o ${meta.id}.gene_abundance.tsv \\
        ${bam}
    """
}
