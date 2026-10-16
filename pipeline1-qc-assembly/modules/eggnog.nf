// Annotates the catalog representatives once. The catalog is nucleotide
// (see MMSEQS2), so --itype CDS makes eggNOG-mapper translate the genes;
// annotation rows keep the catalog gene IDs used in the abundance tables.
process EGGNOG {

    tag "${genes.baseName}"

    label 'process_high'

    // publishDir "${params.outdir}/10_gene_catalog/eggnog", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/eggnogmapper/2.1.12/eggnogmapper_2.1.12.sif"

    containerOptions "--bind ${params.eggnog_db}:/data/eggnog:ro"

    input:
    path(genes)

    output:
    path("${genes.baseName}.emapper.annotations"), emit: annotations
    path("${genes.baseName}.emapper.seed_orthologs"), emit: orthologs

    script:
    """
    emapper.py \\
        -i ${genes} \\
        --itype CDS \\
        --translate \\
        -o ${genes.baseName} \\
        --data_dir /data/eggnog \\
        --cpu ${task.cpus} \\
        -m diamond \\
        --tax_scope auto \\
        --go_evidence non-electronic \\
        --target_orthologs all \\
        --seed_ortholog_evalue 0.001 \\
        --seed_ortholog_score 60 \\
        --query_cover 20 \\
        --subject_cover 0
    """
}

// Joins the per-chunk eggNOG outputs: one column header, then all rows.
// "##" comment lines (run information) are dropped.
process EGGNOG_MERGE {

    tag "gene_catalog"

    label 'process_low'

    publishDir "${params.outdir}/10_gene_catalog/eggnog", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/eggnogmapper/2.1.12/eggnogmapper_2.1.12.sif"

    input:
    path(annotations)
    path(orthologs)

    output:
    path("gene_catalog.emapper.annotations"), emit: annotations
    path("gene_catalog.emapper.seed_orthologs"), emit: orthologs

    script:
    """
    awk '/^#[^#]/ && !header_done { print; header_done = 1; next } !/^#/ { print }' \\
        ${annotations.join(' ')} > gene_catalog.emapper.annotations

    awk '/^#[^#]/ && !header_done { print; header_done = 1; next } !/^#/ { print }' \\
        ${orthologs.join(' ')} > gene_catalog.emapper.seed_orthologs
    """
}