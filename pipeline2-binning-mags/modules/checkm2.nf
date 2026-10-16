process CHECKM2_FINAL {
    // Final quality report on the dereplicated MAG catalog. Per-sample
    // quality (used by dRep) comes from Binette, which runs CheckM2 itself.
    label 'process_high'
    publishDir "${params.outdir}/04_checkm2_final", mode: 'copy'

    // container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/checkm2/0.0.1/checkm2_0.0.1.sif"   
    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/binette/1.2.1/binette_v1.2.1.sif" // CheckM2 1.0.1

    input:
    path bins_dir
    path checkm2_db

    output:
    path "final_mags_checkm2.tsv", emit: results

    script:
    """
    checkm2 predict \\
        --input ${bins_dir} \\
        --output-directory checkm2_out \\
        --threads ${task.cpus} \\
        -x fa \\
        --database_path ${checkm2_db}

    cp checkm2_out/quality_report.tsv final_mags_checkm2.tsv
    """
}