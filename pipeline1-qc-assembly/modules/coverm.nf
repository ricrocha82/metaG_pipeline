process COVERM_CONTIG {
    // Replaces the old minimap2+samtools step: CoverM handles read mapping
    // internally (via minimap2 under the hood) and reports contig-level
    // coverage/abundance directly. We ask it to retain the BAM file too,
    // since the gene abundance step (featureCounts) needs it.
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/08_coverm/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/coverm/0.7.0/coverm_0.7.0_miniconda24.sif"

    input:
    tuple val(meta), path(assembly), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}.contig_abundance.tsv"), emit: abundance
    tuple val(meta), path("bam_cache/*.bam"), emit: bam

    script:
    """
    mkdir -p bam_cache
    coverm contig \\
        --coupled ${reads1} ${reads2} \\
        --reference ${assembly} \\
        --methods count mean length covered_fraction trimmed_mean rpkm tpm \\
        --bam-file-cache-directory bam_cache \\
        --threads ${task.cpus} \\
        --output-file ${meta.id}.contig_abundance.tsv
    """
}



// Maps each sample's clean reads to the nucleotide gene catalog.
// --min-read-percent-identity 95 matches the catalog's 95% clustering, so a
// read is only counted for a gene it could belong to.
process COVERM_GENE_CATALOG {

    tag "$meta.id"

    label 'process_high'

    publishDir { "${params.outdir}/10_gene_catalog/abundance/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/coverm/0.7.0/coverm_0.7.0_miniconda24.sif"

    input:
    tuple val(meta), path(reads1), path(reads2), path(gene_catalog)

    output:
    tuple val(meta), path("${meta.id}.gene_abundance.tsv"), emit: abundance

    script:
    """
    coverm contig \\
        --coupled ${reads1} ${reads2} \\
        --reference ${gene_catalog} \\
        --methods count mean length covered_fraction rpkm tpm \\
        --min-read-percent-identity 95 \\
        --threads ${task.cpus} \\
        --output-file ${meta.id}.gene_abundance.tsv
    """
}