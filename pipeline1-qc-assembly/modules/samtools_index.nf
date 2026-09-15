process SAMTOOLS_INDEX {
    // CoverM's cached BAM is coordinate-sorted internally, but has no .bai
    // yet -- needed by featureCounts/downstream tools.
    tag "$meta.id"
    label 'process_low'

    container "/fs/project/PAS1117/modules/singularity/SAMtools-1.23.1.sif"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("*.bai"), emit: bam_bai

    script:
    """
    samtools index ${bam}
    """
}
