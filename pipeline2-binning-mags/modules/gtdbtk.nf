process GTDBTK {
    // GTDB-Tk 2.6.1 with GTDB release 226. From 2.4 the ANI screen uses the
    // skani sketches in the reference package, so no --mash_db is needed.
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/05_gtdbtk/${meta.id}" }, mode: 'copy'

    container "/fs/ess/PAS0822/modules/singularity/gtdbtk_2.6.1.sif"

    input:
    tuple val(meta), path(bins_dir)
    path gtdbtk_db

    output:
    tuple val(meta), path("${meta.id}_gtdbtk"), emit: results
    path "${meta.id}_gtdbtk/gtdbtk.*.summary.tsv", emit: summaries

    script:
    """
    export GTDBTK_DATA_PATH=\$(readlink -f ${gtdbtk_db})
    gtdbtk classify_wf \\
        --genome_dir ${bins_dir} \\
        --out_dir ${meta.id}_gtdbtk \\
        --extension fa \\
        --cpus ${task.cpus} \\
        --force
    """
}
