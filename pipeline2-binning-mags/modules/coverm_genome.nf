
process COVERM_GENOME_MATRIX {
    // Replaces the earlier per-sample COVERM_GENOME: now maps ALL
    // samples' reads against the dereplicated representative genome set
    // in one go, producing a proper MAG x sample abundance matrix --
    // which is what dereplicated MAGs actually need (their abundance
    // profile across every sample, not just the one they happened to be
    // binned from).
    //
    // R1/R2 pairing comes directly from the input channel tuple (same
    // pairing used for assembly and COVERM_MAPPING earlier in the
    // pipeline) rather than being re-derived from filenames via ls/glob
    // -- avoids the earlier bug where raw-read naming (*_R1_001.fastq.gz)
    // didn't match the assumed *_clean_1.fastq.gz / *_1.fastq.gz pattern.
    label 'process_high'
    publishDir "${params.binning_outdir}/05_mag_abundance", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/CoverM-0.7.0.sif"

    input:
    path representatives
    val coupled_reads

    output:
    path "mag_abundance_matrix.tsv", emit: abundance

    script:
    """
    mkdir -p bam_cache

    echo "Paired samples for CoverM (verify these look correct):"
    echo "${coupled_reads}"

    coverm genome \\
        --coupled ${coupled_reads} \\
        --genome-fasta-directory ${representatives} \\
        -x fa \\
        --methods relative_abundance mean trimmed_mean \\
        --threads ${task.cpus} \\
        --min-read-percent-identity .95 --min-read-aligned-percent .75 --discard-unmapped \\
        --bam-file-cache-directory bam_cache \\
        --output-file mag_abundance_matrix.tsv
    """
}