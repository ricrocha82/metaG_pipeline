// Non-redundant gene catalog across all samples (Maranga et al. 2023 clustered
// complete genes at 95% identity with CD-HIT). Clusters the Prodigal *nucleotide* genes,
// so the same catalog can be used for read mapping (DNA) and for eggNOG-mapper
// (--itype CDS translates them).
// --cluster-mode 2 is greedy like CD-HIT (the longest sequence becomes the
// representative); --cov-mode 1 -c 0.90 requires the representative to cover
// 90% of each member, like CD-HIT's -aS 0.9.
// Only complete genes (Prodigal partial=00: start and stop codon present) are
// clustered, as in the paper; fragments at contig ends are left out.
process MMSEQS2 {

    tag "gene_catalog"

    label 'process_high'

    publishDir "${params.outdir}/10_gene_catalog/mmseqs2", mode: 'copy'

    // eggNOG-mapper image (includes MMseqs2 18): the CoMS mmseqs2 image lacks
    // `ps`, which Nextflow needs inside containers to collect task metrics.
    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/eggnogmapper/2.1.12/eggnogmapper_2.1.12.sif"

    input:
    path(genes)

    output:
    path("gene_catalog_rep_seq.fasta"), emit: representatives
    path("gene_catalog_cluster.tsv"), emit: clusters

    script:
    """
    cat ${genes.join(' ')} \\
        | awk '/^>/ { keep = (\$0 ~ /;partial=00;/) } keep' \\
        > gene_catalog_complete.fna

    mmseqs easy-cluster \\
        gene_catalog_complete.fna \\
        gene_catalog \\
        tmp \\
        --min-seq-id 0.95 \\
        -c 0.90 \\
        --cov-mode 1 \\
        --cluster-mode 2 \\
        --threads ${task.cpus}

    rm -rf tmp gene_catalog_complete.fna
    """
}
