process VIRSORTER2 {
    // --keep-original-seq: report whole contigs and let CheckV do the
    // host-region trimming, as in the CoMS SOP.
    // No --db-dir: the CoMS VirSorter2 image ships with its database.
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/02_virsorter2/${meta.id}" }, mode: 'copy',
        pattern: "vs2/final-viral-*", saveAs: { filename -> file(filename).name }

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/VirSorter2/2.2.3/VirSorter2-2.2.3.sif"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("vs2/final-viral-combined.fa"), emit: virus
    tuple val(meta), path("vs2/final-viral-score.tsv"),   emit: scores
    tuple val(meta), path("vs2/final-viral-boundary.tsv"), emit: boundary

    script:
    """
    virsorter run \\
        --keep-original-seq \\
        -i ${assembly} \\
        -w vs2 \\
        --include-groups ${params.vs2_groups} \\
        --min-length ${params.min_contig_length} \\
        --min-score ${params.vs2_min_score} \\
        -j ${task.cpus} \\
        all
    """
}
