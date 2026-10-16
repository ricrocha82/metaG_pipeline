process DREP {
    // Cross-sample dereplication: collapses genomes found independently in
    // several samples into one representative each (secondary ANI
    // params.drep_secondary_ani), chosen by CheckM2 completeness and
    // contamination from genomeInfo. Genomes below params.drep_min_completeness
    // or above params.drep_max_contamination are removed first.
    //
    // --ignoreGenomeQuality (params.ignore_genome_quality, default false):
    // dRep assumes fractions instead of percentages when every contamination
    // value is < 1 and refuses to run. Only happens on tiny test runs; it
    // also disables the quality filter and quality-based selection.
    label 'process_high'
    publishDir "${params.outdir}/03_drep", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/drep/3.4.0/drep_3.4.0.sif"

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
        -comp ${params.drep_min_completeness} \\
        -con ${params.drep_max_contamination} \\
        -pa 0.9 -sa ${params.drep_secondary_ani} -nc 0.1 \\
        -cm larger \\
        -p ${task.cpus}
    """
}
