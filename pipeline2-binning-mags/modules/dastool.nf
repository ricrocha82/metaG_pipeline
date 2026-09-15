process FASTA_TO_CONTIG2BIN {
    tag "${meta.id}:${binner_label}"
    label 'process_low'

    input:
    tuple val(meta), path(bins_dir), val(binner_label)

    output:
    tuple val(meta), path("${binner_label}_contig2bin.tsv"), emit: table

    script:
    """
    Fasta_to_Contig2Bin.sh -i ${bins_dir} -e fa > ${binner_label}_contig2bin.tsv
    """
}

process DASTOOL {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/02_dastool/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/DAS_Tool-1.1.7.sif"

    input:
    tuple val(meta), path(assembly), path(metabat2_tsv), path(comebin_tsv), path(semibin2_tsv)

    output:
    tuple val(meta), path("${meta.id}_dastool_bins"), emit: bins
    tuple val(meta), path("${meta.id}_DASTool_summary.tsv"), emit: summary

    script:
    """
    DAS_Tool \\
        -i ${metabat2_tsv},${comebin_tsv},${semibin2_tsv} \\
        -l metabat2,comebin,semibin2 \\
        -c ${assembly} \\
        -o ${meta.id} \\
        -t ${task.cpus} \\
        --write_bins \\
        --score_threshold 0.5

    mkdir -p ${meta.id}_dastool_bins
    cp ${meta.id}_DASTool_bins/*.fa ${meta.id}_dastool_bins/ 2>/dev/null || true
    cp ${meta.id}_DASTool_summary.tsv . 2>/dev/null || true
    """
}