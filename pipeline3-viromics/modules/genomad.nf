process GENOMAD {
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/02_genomad/${meta.id}" }, mode: 'copy',
        pattern: "genomad/*_summary/*", saveAs: { filename -> file(filename).name }

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/genomad/1.8.0/genomad_1.8.0.sif"

    input:
    tuple val(meta), path(assembly)
    path db

    output:
    // meta.id matches the file name geNomad derives from the input (${meta.id}.fna)
    tuple val(meta), path("genomad/${meta.id}_summary/${meta.id}_virus.fna"),                                  emit: virus
    // Scores for every contig (not only those passing --min-score); needed
    // for the VS2-exclusive rule in --detect both.
    tuple val(meta), path("genomad/${meta.id}_aggregated_classification/${meta.id}_aggregated_classification.tsv"), emit: scores
    tuple val(meta), path("genomad/${meta.id}_summary/${meta.id}_virus_summary.tsv"),                          emit: summary

    script:
    """
    # geNomad names its outputs after the input file; link it as <sample>.fna
    ln -s ${assembly} ${meta.id}.fna

    genomad end-to-end \\
        --cleanup \\
        --threads ${task.cpus} \\
        --min-score ${params.genomad_min_score} \\
        ${meta.id}.fna \\
        genomad \\
        ${db}
    """
}
