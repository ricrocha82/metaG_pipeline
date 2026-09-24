process GENE_ABUNDANCE {
    tag "$meta.id"
    label 'process_medium'
    publishDir { "${params.outdir}/10_gene_abundance/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/feature-counts/2.0.0/feature-counts_2.0.0.sif"

    input:
    tuple val(meta), path(gff), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.gene_abundance.tsv"), emit: abundance

    script:
    // featureCounts 2.0.0: -p alone counts read pairs (fragments).
    // --countReadPairs only exists from 2.0.2, where -p no longer implies it.
    """
    featureCounts -T ${task.cpus} -p \\
        -t CDS -g ID \\
        -a ${gff} -F GFF \\
        -o ${meta.id}.gene_abundance.tsv \\
        ${bam}
    """
}
