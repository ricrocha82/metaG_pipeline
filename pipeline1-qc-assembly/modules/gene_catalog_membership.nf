// One row per original Prodigal gene: its catalog representative, contig,
// source sample and coordinates. Links the catalog back to contigs, which the
// pipeline2 master table uses to connect genes to MAGs.
// build_gene_membership.py is in bin/, which Nextflow puts on PATH.
process GENE_CATALOG_MEMBERSHIP {

    tag "gene_catalog"

    label 'process_low'

    publishDir "${params.outdir}/10_gene_catalog/membership", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/eggnogmapper/2.1.12/eggnogmapper_2.1.12.sif"

    input:
    path(clusters)
    path(gff_files)

    output:
    path("gene_catalog_membership.tsv"), emit: membership

    script:
    """
    build_gene_membership.py \\
        --clusters ${clusters} \\
        --gff ${gff_files.join(' ')} \\
        --output gene_catalog_membership.tsv
    """
}
