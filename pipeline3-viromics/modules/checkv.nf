process CHECKV {
    // Run once on the pooled candidates from all samples.
    label 'process_high'
    publishDir "${params.outdir}/04_checkv", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/CheckV/0.8.1/CheckV-0.8.1.sif"

    input:
    path candidates
    path db

    output:
    path "checkv/quality_summary.tsv", emit: quality_summary
    path "checkv/viruses.fna",         emit: viruses
    path "checkv/proviruses.fna",      emit: proviruses
    path "checkv/*.tsv",               emit: tables

    script:
    """
    checkv end_to_end \\
        ${candidates} \\
        checkv \\
        -d ${db} \\
        -t ${task.cpus}
    """
}

process FILTER_VIRAL_CONTIGS {
    label 'process_low'
    publishDir "${params.outdir}/05_viral_contigs", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/CheckV/0.8.1/CheckV-0.8.1.sif"

    input:
    path quality_summary
    path viruses
    path proviruses
    path candidates_table

    output:
    path "viral_contigs.fna",         emit: fasta
    path "viral_contigs.tsv",         emit: table
    // Viral sequences that failed a rule (e.g. < min_contig_length), with the reason
    path "removed_viral_contigs.fna", emit: removed_fasta
    path "removed_viral_contigs.tsv", emit: removed_table

    script:
    """
    filter_viral_contigs.py \\
        --quality-summary ${quality_summary} \\
        --viruses ${viruses} \\
        --proviruses ${proviruses} \\
        --candidates ${candidates_table} \\
        --mode ${params.detect} \\
        --min-length ${params.min_contig_length} \\
        --exclusive-min-length ${params.exclusive_min_length} \\
        --exclusive-min-viral-genes ${params.exclusive_min_viral_genes} \\
        --exclusive-min-genomad-score ${params.exclusive_min_genomad_score} \\
        --out-prefix viral_contigs
    """
}
