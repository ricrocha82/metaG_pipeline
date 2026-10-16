process GTDBTK_DE_NOVO {
    // iPHoP add_to_db needs GTDB-Tk de novo trees for bacteria AND archaea
    // in the same output folder, so both runs live in one task.
    // The archaeal run fails when there are no archaeal MAGs; that is
    // expected and tolerated (as in the CoMS SOP).
    label 'process_high'
    publishDir "${params.outdir}/09_iphop", mode: 'copy'

    container params.gtdbtk_sif
    // containerOptions "--bind ${params.gtdbtk_db}:/refdata"

    input:
    path mags, stageAs: 'mags/*'

    output:
    path "gtdbtk_de_novo", emit: gtdbtk

    script:
    """
    export GTDBTK_DATA_PATH=\$(readlink -f ${params.gtdbtk_db})
    
    gtdbtk de_novo_wf \\
        --genome_dir mags \\
        --bacteria \\
        --outgroup_taxon p__Patescibacteriota \\
        --out_dir gtdbtk_de_novo \\
        --cpus ${task.cpus} \\
        --force \\
        -x ${params.host_mags_ext}

    if ! gtdbtk de_novo_wf \\
        --genome_dir mags \\
        --archaea \\
        --outgroup_taxon p__Altiarchaeota \\
        --out_dir gtdbtk_de_novo \\
        --cpus ${task.cpus} \\
        --force \\
        -x ${params.host_mags_ext}
    then
        echo "GTDB-Tk archaeal de novo run failed (usually: no archaeal MAGs). Continuing with bacteria only." >&2
    fi
    """
}

process IPHOP_ADD_TO_DB {
    // The custom DB is a full iPHoP DB (default hosts + your MAGs), so
    // predictions are run against it alone. storeDir keeps it outside
    // work/ and skips the rebuild on later runs; pass it with --iphop_db
    // to reuse it without --host_mags.
    label 'process_high'
    storeDir "${params.outdir}/09_iphop/custom_db"

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/iphop/1.4.1/iphop_1.4.1.sif"
    containerOptions '--cleanenv'

    input:
    path mags, stageAs: 'mags/*'
    path gtdbtk
    path default_db

    output:
    path "iphop_custom_db", emit: db

    script:
    """
    # The CoMS iPHoP image has no iphop on PATH: its runscript sets up the
    # environment (the SOP's `apptainer run ... add_to_db`), so call it directly.
    /.singularity.d/runscript add_to_db \\
        --fna_dir mags \\
        --gtdb_dir ${gtdbtk} \\
        --out_dir iphop_custom_db \\
        --db_dir ${default_db}
    """
}

process IPHOP_PREDICT {
    tag "${chunk.baseName}"
    label 'process_high'
    publishDir { "${params.outdir}/09_iphop/chunks/${chunk.baseName}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/iphop/1.4.1/iphop_1.4.1.sif"
    containerOptions '--cleanenv'

    input:
    path chunk
    path db

    output:
    path "iphop_out/Host_prediction_to_genus_m90.csv",  emit: genus
    path "iphop_out/Host_prediction_to_genome_m90.csv", emit: genome
    path "iphop_out/Detailed_output_by_tool.csv",       emit: detailed

    script:
    """
    # Runscript = the SOP's `apptainer run ... predict` (iphop is not on PATH)
    /.singularity.d/runscript predict \\
        --fa_file ${chunk} \\
        --db_dir ${db} \\
        --out_dir iphop_out \\
        -t ${task.cpus}
    """
}
