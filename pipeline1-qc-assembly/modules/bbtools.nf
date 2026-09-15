// Adapter + quality trimming in one BBDuk pass (kmer adapter trim + quality
// trim), mirroring ricardo's existing bbduk workflow but combined into a
// single call rather than three sequential ones.
process BBDUK_TRIM {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.qc_outdir}/02_bbduk_trim/${meta.id}", mode: 'copy', pattern: "*.stats"

    container "/fs/project/PAS1117/modules/singularity/BBTools-39.31.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path adapters_fasta

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

// Host read removal via bbmap.sh alignment-based filtering, matching
// ricardo's existing script's approach (step 4 in his bbduk.sh) but as its
// own explicit pipeline step.
process BBMAP_HOST_REMOVAL {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.qc_outdir}/03_host_removal/${meta.id}", mode: 'copy', pattern: "*.stats"

    container "/fs/project/PAS1117/modules/singularity/BBTools-39.31.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path host_fasta

    output:
    tuple val(meta), path("${meta.id}_host_removed_1.fastq.gz"), path("${meta.id}_host_removed_2.fastq.gz"), emit: reads
    path "${meta.id}_hostRemoval.stats", emit: stats

    script:
    """
    bbmap.sh -Xmx${task.memory.toGiga()}g threads=${task.cpus} \\
        minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \\
        overwrite=t \\
        ref=${host_fasta} \\
        in1=${reads1} in2=${reads2} \\
        outu1=${meta.id}_host_removed_1.fastq.gz outu2=${meta.id}_host_removed_2.fastq.gz \\
        outm1=${meta.id}_host_matched_1.fastq.gz outm2=${meta.id}_host_matched_2.fastq.gz \\
        statsfile=${meta.id}_hostRemoval.stats \\
        nodisk
    """
}

// PhiX removal -- a *dedicated* PhiX reference, not the adapters.fa
// (ricardo's original script pointed ref= at adapters.fa for this step,
// which looks like a copy-paste leftover -- fixed here).
process BBDUK_PHIX_REMOVAL {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.qc_outdir}/03b_phix_removal/${meta.id}", mode: 'copy', pattern: "*.stats"
    publishDir "${params.qc_outdir}/04a_clean_reads", mode: 'copy', pattern: "*_clean_{1,2}.fastq.gz"

    container "/fs/project/PAS1117/modules/singularity/BBTools-39.31.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path phix_fasta

    output:
    tuple val(meta), path("${meta.id}_clean_1.fastq.gz"), path("${meta.id}_clean_2.fastq.gz"), emit: reads
    path "${meta.id}_PhiXFiltering.stats", emit: stats

    script:
    """
    bbduk.sh -Xmx${task.memory.toGiga()}g threads=${task.cpus} \\
        overwrite=t \\
        in1=${reads1} in2=${reads2} \\
        out1=${meta.id}_clean_1.fastq.gz out2=${meta.id}_clean_2.fastq.gz \\
        ref=${phix_fasta} \\
        k=31 hdist=1 \\
        stats=${meta.id}_PhiXFiltering.stats statscolumns=5
    """
}
