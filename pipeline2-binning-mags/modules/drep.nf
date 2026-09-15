process DREP {
    // Cross-sample dereplication: collapses near-identical genomes
    // (default dRep threshold: 99% ANI) found independently across your
    // samples down to one representative each, chosen by genome quality
    // (completeness/contamination from CheckM2).
    //
    // --ignoreGenomeQuality (toggled via params.ignore_genome_quality,
    // default false): dRep has a sanity-check heuristic that misfires
    // with very few genomes -- if EVERY contamination value in the
    // dataset happens to be under 1 (even if genuinely valid, e.g. a
    // real 0.03% contamination), it assumes you supplied fractions
    // instead of percentages and refuses to proceed. This is a known
    // edge case with small test runs (e.g. a single sample/genome), not
    // a real data problem. Set true only for testing; leave false for
    // real multi-sample runs, where CheckM2-based quality filtering is
    // the actual point of supplying genomeInfo at all.
    label 'process_high'
    publishDir "${params.binning_outdir}/03b_drep", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/dRep-3.6.2.sif"

    input:
    path pooled_bins
    path genome_info

    output:
    path "drep_out/dereplicated_genomes", emit: representatives
    path "drep_out/data_tables", emit: tables

    script:
    def quality_flag = params.ignore_genome_quality ? '--ignoreGenomeQuality' : "--genomeInfo ${genome_info}"
    """
    dRep dereplicate drep_out \\
        -g ${pooled_bins}/*.fa \\
        ${quality_flag} \\
        -comp ${params.drep_min_completeness ?: 75} \\
        -con ${params.drep_max_contamination ?: 25} \\
        -pa 0.9 -sa 0.95 -nc 0.1 \\
        -cm larger \\
        -p ${task.cpus}
    """
}
