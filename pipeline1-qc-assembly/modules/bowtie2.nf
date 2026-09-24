// Builds a Bowtie2 index once. storeDir keeps it outside the work directory,
// keyed by the reference file name, so later runs (or a new outdir pointing to
// the same bowtie2_index_dir) reuse it instead of rebuilding the host genome.
process BOWTIE2_BUILD {
    tag "$index_name"
    label 'process_high'
    storeDir "${params.bowtie2_index_dir ?: "${params.outdir}/bowtie2_indexes"}/${reference.name}"

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/bowtie2/2.5.1/bowtie2_2.5.1.sif"

    input:
    tuple val(index_name), path(reference)

    output:
    tuple val(index_name), path("${index_name}.*.bt2*"), emit: index

    script:
    """
    bowtie2-build \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${index_name}
    """
}

// Keeps read pairs that do not align concordantly to the reference
// (--un-conc-gz), matching the CoMS preprocessing workflow. Used twice via
// aliases: host removal, then PhiX removal.
process BOWTIE2_REMOVAL {
    tag "${meta.id} (${index_name})"
    label 'process_high'
    publishDir { "${params.outdir}/03_${index_name}_removal" }, mode: 'copy', pattern: "*.log"
    publishDir "${params.outdir}/04a_clean_reads", mode: 'copy', pattern: "*_phix_removed_*.fastq.gz",
        saveAs: { filename -> filename.replace('_phix_removed_', '_clean_') }

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/bowtie2/2.5.1/bowtie2_2.5.1.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    tuple val(index_name), path(index_files)

    output:
    tuple val(meta), path("${meta.id}_${index_name}_removed_1.fastq.gz"), path("${meta.id}_${index_name}_removed_2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_${index_name}_mapped_1.fastq.gz"), path("${meta.id}_${index_name}_mapped_2.fastq.gz"), emit: mapped_reads
    path "${meta.id}_${index_name}.bowtie2.log", emit: log

    script:
    """
    bowtie2 \\
        --end-to-end --sensitive \\
        -p ${task.cpus} \\
        -x ${index_name} \\
        -1 ${reads1} -2 ${reads2} \\
        --un-conc-gz ${meta.id}_${index_name}_removed_%.fastq.gz \\
        --al-conc-gz ${meta.id}_${index_name}_mapped_%.fastq.gz \\
        -S /dev/null \\
        2> ${meta.id}_${index_name}.bowtie2.log
    """
}
