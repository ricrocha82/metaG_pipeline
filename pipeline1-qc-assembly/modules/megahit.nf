process MEGAHIT {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.qc_outdir}/06_megahit/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/MEGAHIT-1.2.9.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}.contigs.fa"), emit: assembly

    script:
    """
    megahit -1 ${reads1} -2 ${reads2} \\
        -t ${task.cpus} \\
        --min-contig-len ${params.qc_min_contig_size ?: 1000} \\
        -o megahit_out

    cp megahit_out/final.contigs.fa ${meta.id}.contigs.fa
    """
}
