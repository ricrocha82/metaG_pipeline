// Adapter + quality trimming in one BBDuk pass (kmer adapter trim + quality trim).
process BBDUK_TRIM {
    tag "$meta.id"
    label 'process_medium'
    publishDir { "${params.outdir}/02_bbduk_trim/${meta.id}" }, mode: 'copy', pattern: "*.stats"

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/BBTools/38.69/BBTools-38.69.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    // val, not path: either BBDuk's keyword 'adapters' (reference bundled in
    // the container) or a fasta path on a filesystem visible to the container.
    val adapters_fasta

    output:
    tuple val(meta), path("${meta.id}_trimmed_1.fastq.gz"), path("${meta.id}_trimmed_2.fastq.gz"), emit: reads
    path "${meta.id}_adapterTrimming.stats", emit: stats

    script:
    """
    bbduk.sh -Xmx${task.memory.toGiga()}g threads=${task.cpus} \\
        overwrite=t \\
        in1=${reads1} in2=${reads2} \\
        out1=${meta.id}_trimmed_1.fastq.gz out2=${meta.id}_trimmed_2.fastq.gz \\
        ref=${adapters_fasta} \\
        ktrim=r k=23 mink=11 hdist=1 \\
        qtrim=rl trimq=14 minlength=30 maq=20 maxns=0 tpe tbo \\
        refstats=${meta.id}_adapterTrimming.stats statscolumns=5
    """
}
