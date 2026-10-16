process COVERM_GENOME_MATRIX {
    // Maps ALL samples' reads against the dereplicated representative
    // genomes in one CoverM call, producing a MAG x sample abundance matrix.
    // A MAG counts as present in a sample only if at least
    // params.coverm_min_covered_fraction of its length is covered, which
    // limits false presence calls from reads on conserved regions.
    //
    // R1/R2 pairing comes from the samplesheet (same pairing used for
    // COVERM_MAPPING), not re-derived from filenames.
    label 'process_high'
    publishDir "${params.outdir}/06_mag_abundance", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/coverm/0.7.0/coverm_0.7.0_miniconda24.sif"

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
        --methods relative_abundance count mean trimmed_mean \\
        --threads ${task.cpus} \\
        --min-read-percent-identity .95 --min-read-aligned-percent .75 \\
        --min-covered-fraction ${params.coverm_min_covered_fraction} \\
        --discard-unmapped \\
        --bam-file-cache-directory bam_cache \\
        --output-file mag_abundance_matrix.tsv
    """
}
