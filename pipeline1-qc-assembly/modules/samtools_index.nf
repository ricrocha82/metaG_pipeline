process SAMTOOLS_INDEX {
    // CoverM's cached BAM is coordinate-sorted internally, but has no .bai
    // yet -- needed by featureCounts/downstream tools.
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
