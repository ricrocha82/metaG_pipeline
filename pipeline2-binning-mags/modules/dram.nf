// DRAM 1.4.6 (CoMS image without KEGG). As in the CoMS workflow, each MAG
// is annotated in its own job, then annotations are merged and distilled.
// The CoMS DRAM databases are bound at /databases/DRAM, where the image's
// DRAM config expects them.

process DRAM_ANNOTATE {
    tag "$meta.id"
    label 'process_high'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/DRAM/1.4.6_no_kegg/DRAM-1.4.6_nokegg.sif"
    containerOptions "--bind ${params.dram_db}:/databases/DRAM --writable-tmpfs --fakeroot"

    input:
    tuple val(meta), path(mag)
    path gtdbtk_summaries
    path checkm2_report

    output:
    tuple val(meta), path("${meta.id}_dram"), emit: annotations

    script:
    // GTDB-Tk writes one summary per domain (bac120, ar53); DRAM accepts
    // --gtdb_taxonomy once per file.
    def taxonomy_args = gtdbtk_summaries.collect { summary -> "--gtdb_taxonomy ${summary}" }.join(' ')
    """
    DRAM.py annotate \\
        -i ${mag} \\
        -o ${meta.id}_dram \\
        ${taxonomy_args} \\
        --checkm_quality ${checkm2_report} \\
        --threads ${task.cpus}
    """
}

process DRAM_MERGE {
    label 'process_low'
    publishDir "${params.outdir}/07_dram", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/DRAM/1.4.6_no_kegg/DRAM-1.4.6_nokegg.sif"
    containerOptions "--bind ${params.dram_db}:/databases/DRAM --writable-tmpfs --fakeroot"

    input:
    path annotation_dirs

    output:
    path "dram_merged", emit: merged

    script:
    """
    DRAM.py merge_annotations -i '*_dram' -o dram_merged
    """
}

process DRAM_DISTILL {
    label 'process_medium'
    publishDir "${params.outdir}/07_dram", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/DRAM/1.4.6_no_kegg/DRAM-1.4.6_nokegg.sif"
    containerOptions "--bind ${params.dram_db}:/databases/DRAM --writable-tmpfs --fakeroot"

    input:
    path merged

    output:
    path "dram_distill", emit: distill

    script:
    """
    # tRNA/rRNA tables are only written when DRAM found any.
    OPTIONAL_ARGS=""
    if [ -f ${merged}/trnas.tsv ]; then OPTIONAL_ARGS="\$OPTIONAL_ARGS --trna_path ${merged}/trnas.tsv"; fi
    if [ -f ${merged}/rrnas.tsv ]; then OPTIONAL_ARGS="\$OPTIONAL_ARGS --rrna_path ${merged}/rrnas.tsv"; fi

    DRAM.py distill \\
        -i ${merged}/annotations.tsv \\
        -o dram_distill \\
        \$OPTIONAL_ARGS
    """
}
