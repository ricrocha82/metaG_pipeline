process PRODIGAL {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.qc_outdir}/09_prodigal/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/Prodigal-2.6.3.sif"

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

process BAKTA {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.qc_outdir}/09_bakta/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/Bakta-1.12.0.sif"

    input:
    tuple val(meta), path(assembly)
    val bakta_db

    output:
    tuple val(meta), path("${meta.id}_bakta"), emit: results
    tuple val(meta), path("${meta.id}_bakta/${meta.id}.gff3"), emit: gff
    tuple val(meta), path("${meta.id}_bakta/${meta.id}.faa"), emit: proteins
    tuple val(meta), path("${meta.id}_bakta/${meta.id}.ffn"), emit: genes
    tuple val(meta), path("${meta.id}_bakta/${meta.id}.tsv"), emit: annotations
    tuple val(meta), path("${meta.id}_bakta/${meta.id}.json"), emit: json

    script:
    """
    test -d "${bakta_db}" || {
        echo "ERROR: Bakta database not found: ${bakta_db}" >&2
        exit 1
    }

    bakta \
        --db "${bakta_db}" \
        --output "${meta.id}_bakta" \
        --prefix "${meta.id}" \
        --threads ${task.cpus} \
        --force \
        "${assembly}"
    """
}