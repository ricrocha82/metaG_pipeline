process MEGAHIT {
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/06_megahit/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/MEGAHIT/1.2.8/MEGAHIT-1.2.8.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: assembly

    script:
    // MEGAHIT names contigs k141_N in every sample, so IDs collide when
    // contigs from several samples are pooled downstream. Prefix the sample
    // ID and drop the description (flag=, multi=, len=) after the first space.
    """
    megahit -1 ${reads1} -2 ${reads2} \\
        -t ${task.cpus} \\
        --min-contig-len ${params.min_contig_size ?: 1000} \\
        -o megahit_out

    awk -v sample="${meta.id}" '/^>/ { sub(/ .*/, ""); sub(/^>/, ">" sample "_") } { print }' \\
        megahit_out/final.contigs.fa > ${meta.id}.contigs.fa
    """
}
