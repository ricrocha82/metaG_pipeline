process PRODIGAL {
    tag "$meta.id"
    label 'process_low'

    publishDir { "${params.outdir}/09_prodigal/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/prodigal/2.6.3/prodigal_2.6.3.sif"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.genes.gff"), emit: gff
    tuple val(meta), path("${meta.id}.genes.fna"), emit: nt_fasta
    tuple val(meta), path("${meta.id}.genes.faa"), emit: aa_fasta

    script:
    """
    prodigal \
        -i ${assembly} \
        -p meta \
        -f gff \
        -o ${meta.id}.genes.gff \
        -d ${meta.id}.genes.fna \
        -a ${meta.id}.genes.faa
    """
}

process PROKKA {
    tag "$meta.id"
    label 'process_medium'

    publishDir { "${params.outdir}/09_prokka/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/prokka/1.14.6/prokka_1.14.6.sif"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_prokka"), emit: results
    tuple val(meta), path("${meta.id}_prokka/${meta.id}.gff"), emit: gff
    tuple val(meta), path("${meta.id}_prokka/${meta.id}.faa"), emit: proteins
    tuple val(meta), path("${meta.id}_prokka/${meta.id}.ffn"), emit: genes
    tuple val(meta), path("${meta.id}_prokka/${meta.id}.tsv"), emit: annotations

    script:
    // Prokka rejects contig IDs longer than 37 characters; MEGAHIT contigs
    // are named <sample>_k141_N, so keep sample IDs short.
    """
    prokka \\
        --metagenome \\
        --outdir ${meta.id}_prokka \\
        --prefix ${meta.id} \\
        --locustag ${meta.id} \\
        --cpus ${task.cpus} \\
        --force \\
        ${assembly}
    """
}
