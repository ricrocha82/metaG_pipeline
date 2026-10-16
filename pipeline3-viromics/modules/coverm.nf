process COVERM_VOTU {
    // Read mapping with minimap2 (through CoverM) against the vOTU set.
    // Thresholds follow the CoMS SOP: reads >= 95% identity over >= 75%
    // of their length; a vOTU with < 70% of its length covered is set to 0.
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/07_coverm/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/coverm/0.7.0/coverm_0.7.0_miniconda24.sif"

    input:
    tuple val(meta), path(reads1), path(reads2), path(votus)

    output:
    tuple val(meta), path("${meta.id}.coverm.tsv"),                 emit: abundance
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam

    script:
    """
    coverm contig \\
        --coupled ${reads1} ${reads2} \\
        --reference ${votus} \\
        --mapper minimap2-sr \\
        --min-read-percent-identity ${params.coverm_min_identity} \\
        --min-read-aligned-percent ${params.coverm_min_aligned} \\
        --min-covered-fraction 0 \\
        -m mean trimmed_mean covered_bases covered_fraction variance length count reads_per_base rpkm tpm \\
        --output-format sparse \\
        --bam-file-cache-directory bam_cache \\
        --threads ${task.cpus} \\
        --output-file ${meta.id}.coverm.tsv

    mv bam_cache/*.bam ${meta.id}.bam
    samtools index ${meta.id}.bam
    """
}

process VOTU_ABUNDANCE_TABLES {
    label 'process_low'
    publishDir "${params.outdir}/08_abundance", mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/CheckV/0.8.1/CheckV-0.8.1.sif"

    input:
    path coverm_tables

    output:
    path "vOTU_*.tsv", emit: tables

    script:
    """
    votu_abundance_tables.py \\
        --min-covered-fraction ${params.coverm_min_covered_fraction} \\
        --out-prefix vOTU \\
        ${coverm_tables}
    """
}
