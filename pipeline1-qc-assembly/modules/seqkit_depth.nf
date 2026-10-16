// Sequencing depth of the clean reads (after trimming, host and PhiX removal).
// These are the reads mapped by CoverM downstream, so their total bases are
// the right denominator for depth-normalizing coverage (e.g. trimmed mean / Gbp).
process SEQKIT_DEPTH {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}/04b_read_depth", mode: 'copy', pattern: "*.seqkit_stats.tsv"

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/seqkit/2.7.0/seqkit_2.7.0.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}.depth.tsv"), emit: depth
    path "${meta.id}.seqkit_stats.tsv", emit: stats

    script:
    """
    seqkit stats -T -j ${task.cpus} ${reads1} ${reads2} > ${meta.id}.seqkit_stats.tsv

    # Sum R1 + R2 (columns: file format type num_seqs sum_len ...).
    # %.0f instead of %d avoids integer overflow in mawk for large libraries.
    awk -v sample=${meta.id} 'NR > 1 { reads += \$4; bases += \$5 }
        END { printf "%s\\t%.0f\\t%.0f\\t%.4f\\n", sample, reads, bases, bases / 1e9 }' \\
        ${meta.id}.seqkit_stats.tsv > ${meta.id}.depth.tsv
    """
}
