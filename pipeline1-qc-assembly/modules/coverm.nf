process COVERM_CONTIG {
    // Replaces the old minimap2+samtools step: CoverM handles read mapping
    // internally (via minimap2 under the hood) and reports contig-level
    // coverage/abundance directly. We ask it to retain the BAM file too,
    // since the gene abundance step (featureCounts) needs it.
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/08_coverm/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/CoverM-0.7.0.sif"

    input:
    tuple val(meta), path(assembly), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}.contig_abundance.tsv"), emit: abundance
    tuple val(meta), path("bam_cache/*.bam"), emit: bam

    script:
    """
    mkdir -p bam_cache
    coverm contig \\
        --coupled ${reads1} ${reads2} \\
        --reference ${assembly} \\
        --methods mean trimmed_mean covered_fraction \\
        --bam-file-cache-directory bam_cache \\
        --threads ${task.cpus} \\
        --output-file ${meta.id}.contig_abundance.tsv
    """
}
