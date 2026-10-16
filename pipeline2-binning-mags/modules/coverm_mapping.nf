process COVERM_MAPPING {
    // Own mapping step (independent of pipeline1) so this pipeline is
    // self-contained -- just needs assembly + reads. Produces the BAM used
    // for depth calculation across all binners.
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/00_coverm_mapping/${meta.id}" }, mode: 'copy', pattern: "*.tsv"

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/coverm/0.7.0/coverm_0.7.0_miniconda24.sif"

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

process SAMTOOLS_INDEX {
    // CoverM's cached BAM is already coordinate-sorted; it only needs an index.
    tag "$meta.id"
    label 'process_low'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/samtools/1.17/samtools_1.17-2023-06.sif"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam_bai

    script:
    """
    samtools index ${bam}
    """
}
