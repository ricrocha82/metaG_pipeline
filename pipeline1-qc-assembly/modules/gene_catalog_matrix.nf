// Merges the per-sample CoverM tables into gene x sample matrices.
// merge_gene_abundance.py is in bin/, which Nextflow puts on PATH.
process GENE_CATALOG_MATRIX {

    tag "gene_catalog"

    label 'process_low'

    publishDir "${params.outdir}/10_gene_catalog/matrices", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/eggnogmapper/2.1.12/eggnogmapper_2.1.12.sif"

    input:
    path(gene_catalog)
    path(abundance_tables)

    output:
    path("gene_catalog_counts.tsv"), emit: counts
    path("gene_catalog_cpm.tsv"), emit: cpm
    path("gene_catalog_mean_coverage.tsv"), emit: mean_coverage
    path("gene_catalog_tpm.tsv"), emit: tpm

    script:
    """
    merge_gene_abundance.py \\
        --catalog ${gene_catalog} \\
        --tables ${abundance_tables.join(' ')} \\
        --outdir .
    """
}