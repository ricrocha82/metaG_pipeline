// process GTDBTK_PREPARE_DB {
//     // Your provided gtdbtk_db path is a .tar.gz archive; GTDB-Tk needs it
//     // extracted, with GTDBTK_DATA_PATH pointing at the resulting directory.
//     // Cached via storeDir so this only runs once, not on every invocation.
//     label 'process_low'
//     storeDir "${params.binning_outdir}/gtdbtk_db_extracted"

//     input:
//     path db_archive

//     output:
//     path "gtdbtk_data", emit: db_dir

//     script:
//     """
//     mkdir -p gtdbtk_data
//     tar xvzf ${db_archive} -C gtdbtk_data --strip-components=1
//     """
// }

process GTDBTK {
    // CAUTION (flagged before building): gtdbtk_r95_data.tar.gz is from
    // GTDB release 95 (~2020), but GTDB-Tk-2.1.1-PAS1117.sif is a v2.x
    // binary, which expects reference package release ~207+. This will
    // likely fail GTDB-Tk's internal version compatibility check. Verify
    // with: `find /fs/project/PAS1117 -iname "*gtdbtk*"` for a newer
    // release before relying on this step's output.
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.binning_outdir}/04_gtdbtk/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/GTDB-Tk-2.1.1-PAS1117.sif"

    input:
    tuple val(meta), path(bins_dir)
    val gtdbtk_db

    output:
    tuple val(meta), path("${meta.id}_gtdbtk"), emit: results

    script:
    """
    export GTDBTK_DATA_PATH=${gtdbtk_db}
    gtdbtk classify_wf \\
        --genome_dir ${bins_dir} \\
        --out_dir ${meta.id}_gtdbtk \\
        --extension fa \\
        --cpus ${task.cpus} \\
        --force
    """
}
