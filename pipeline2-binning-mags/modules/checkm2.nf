process CHECKM2_DOWNLOAD_DB {
    // Only runs if params.checkm2_db isn't supplied. CheckM2 needs its own
    // DIAMOND reference (uniref100.KO.1.dmnd) -- NOT the CheckM v1
    // CheckMdata format you have at /fs/project/PAS1117/modules/CheckMdata.
    label 'process_low'
    storeDir "${params.binning_outdir}/checkm2_db"

    container "/fs/project/PAS1117/modules/singularity/CheckM2-1.0.1.sif"

    output:
    path "CheckM2_database", emit: db

    script:
    """
    checkm2 database --download --path .
    """
}

process CHECKM2_FINAL {
    // Final quality report on the dereplicated, non-redundant MAG catalog
    // -- these are the numbers that matter for your actual MAG set,
    // versus the earlier per-sample CHECKM2 runs which existed mainly to
    // feed dRep's genomeInfo (quality-based representative selection).
    label 'process_high'
    publishDir "${params.binning_outdir}/06_checkm2_final", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/CheckM2-1.0.1.sif"

    input:
    path bins_dir
    path checkm2_db

    output:
    path "final_mags_checkm2.tsv", emit: results

    script:
    """
    DB_FILE=\$(find -L ${checkm2_db} -name "*.dmnd" | head -n1)
    checkm2 predict \\
        --input ${bins_dir} \\
        --output-directory checkm2_out \\
        --threads ${task.cpus} \\
        -x fa \\
        --database_path \${DB_FILE}

    cp checkm2_out/quality_report.tsv final_mags_checkm2.tsv
    """
}

process CHECKM2 {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.binning_outdir}/03_checkm2/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/CheckM2-1.0.1.sif"

    input:
    tuple val(meta), path(bins_dir)
    path checkm2_db

    output:
    tuple val(meta), path("${meta.id}_checkm2.tsv"), emit: results

    script:
    """
    if [ -n "\$(ls -A ${bins_dir}/*.fa 2>/dev/null)" ]; then
        DB_FILE=\$(find -L ${checkm2_db} -name "*.dmnd" | head -n1)
        checkm2 predict \\
            --input ${bins_dir} \\
            --output-directory checkm2_out \\
            --threads ${task.cpus} \\
            -x fa \\
            --database_path \${DB_FILE}

        cp checkm2_out/quality_report.tsv ${meta.id}_checkm2.tsv
    else
        echo -e "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes" > ${meta.id}_checkm2.tsv
        echo "No bins to score for ${meta.id} -- refinement produced zero bins passing threshold" >&2
    fi
    """
}

process CHECKM2_SING {
    tag "${meta.id}_${meta.binner}"
    label 'process_high'
    publishDir "${params.binning_outdir}/01_binning/${meta.id}/checkm2", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/CheckM2-1.0.1.sif"

    input:
    tuple val(meta), path(bins_dir)
    path checkm2_db

    output:
    tuple val(meta), path("${meta.binner}_checkm2.tsv"), emit: results

    script:
    """
    if [ -n "\$(ls -A ${bins_dir}/*.fa 2>/dev/null)" ]; then
        DB_FILE=\$(find -L ${checkm2_db} -name "*.dmnd" | head -n1)
        checkm2 predict \\
            --input ${bins_dir} \\
            --output-directory checkm2_out \\
            --threads ${task.cpus} \\
            -x fa \\
            --database_path \${DB_FILE}

        cp checkm2_out/quality_report.tsv ${meta.binner}_checkm2.tsv
    else
        echo -e "Name\\tCompleteness\\tContamination\\tCompleteness_Model_Used\\tTranslation_Table_Used\\tCoding_Density\\tContig_N50\\tAverage_Gene_Length\\tGenome_Size\\tGC_Content\\tTotal_Coding_Sequences\\tAdditional_Notes" > ${meta.id}_checkm2.tsv
        echo -e "Name\\t...\\tAdditional_Notes" > ${meta.binner}_checkm2.tsv
        echo "No bins to score for ${meta.id} (${meta.binner}) -- refinement produced zero bins passing threshold" >&2
    fi
    """
}
