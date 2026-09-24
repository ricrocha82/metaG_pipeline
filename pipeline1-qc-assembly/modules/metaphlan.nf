process METAPHLAN {
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/05c_metaphlan/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/metaphlan4/4.2/metaphlan4.2_miniconda24.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path db

    output:
    tuple val(meta), path("${meta.id}.metaphlan_profile.txt"), emit: profile
    path "${meta.id}.metaphlan.mapout.bz2", emit: mapout

    script:
    // --index must name a database already in db: compute nodes have no
    // internet, and without it MetaPhlAn tries to download the latest one.
    // Mates are passed comma-separated; MetaPhlAn profiles them as single reads.
    """
    metaphlan ${reads1},${reads2} \\
        --input_type fastq \\
        -t rel_ab_w_read_stats \\
        --db_dir ${db} \\
        --index ${params.metaphlan_index} \\
        --mapout ${meta.id}.metaphlan.mapout.bz2 \\
        --nproc ${task.cpus} \\
        -o ${meta.id}.metaphlan_profile.txt

    """
}
