process COVERM_MAPPING {
    // Own mapping step (independent of pipeline1) so this pipeline is
    // self-contained -- just needs assembly + reads. Produces the BAM used
    // for depth calculation across all binners.
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/00_coverm_mapping/${meta.id}", mode: 'copy', pattern: "*.tsv"

    container "/fs/project/PAS1117/modules/singularity/CoverM-0.7.0.sif"

    input:
    tuple val(meta), path(assembly), path(reads1), path(reads2)

    output:
    tuple val(meta), path("bam_cache/*.bam"), emit: bam
    tuple val(meta), path("${meta.id}.contig_depth.tsv"), emit: depth

    script:
    """
    mkdir -p bam_cache
    coverm contig \\
        --coupled ${reads1} ${reads2} \\
        --reference ${assembly} \\
        --methods mean \\
        --bam-file-cache-directory bam_cache \\
        --threads ${task.cpus} \\
        --output-file ${meta.id}.contig_depth.tsv
    """
}

process SAMTOOLS_SORT_INDEX {
    tag "$meta.id"
    label 'process_low'

    container "/fs/project/PAS1117/modules/singularity/SAMtools-1.23.1.sif"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam_bai

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam ${bam}
    samtools index ${meta.id}.sorted.bam
    """
}
