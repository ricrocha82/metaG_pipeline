process SINGLEM_PIPE {
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/05b_singlem/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/singlem/100b5/singelm100b5_dmnd_corr.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path metapackage

    output:
    tuple val(meta), path("${meta.id}.otu_table.tsv"), emit: otu_table
    tuple val(meta), path("${meta.id}.profile.tsv"), emit: profile

    script:
    """
    singlem pipe \\
        -1 ${reads1} \\
        -2 ${reads2} \\
        --metapackage ${metapackage} \\
        --otu-table ${meta.id}.otu_table.tsv \\
        -p ${meta.id}.profile.tsv \\
        --threads ${task.cpus}
    """
}

// Fraction of reads that are bacterial or archaeal (the rest being host,
// eukaryotes, phages, etc.). In SingleM 1.0.0beta5 this subcommand is called
// read_fraction; newer SingleM releases renamed it prokaryotic_fraction.
process SINGLEM_READ_FRACTION {
    tag "$meta.id"
    label 'process_low'
    publishDir { "${params.outdir}/05b_singlem/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/singlem/100b5/singelm100b5_dmnd_corr.sif"

    input:
    tuple val(meta), path(reads1), path(reads2), path(profile)
    path metapackage

    output:
    tuple val(meta), path("${meta.id}.read_fraction.tsv"), emit: read_fraction

    script:
    """
    read_fraction_sample=\$(basename ${reads1} .gz)
    awk -F '\\t' -v OFS='\\t' -v sample="\$read_fraction_sample" 'NR > 1 { \$1 = sample } { print }' \\
                        ${profile} > profile_for_read_fraction.tsv

    singlem read_fraction \\
        -1 ${reads1} \\
        -2 ${reads2} \\
        -p profile_for_read_fraction.tsv \\
        --metapackage ${metapackage} \\
        --output-tsv read_fraction.tsv

    awk -F '\\t' -v OFS='\\t' -v sample="${meta.id}" 'NR > 1 { \$1 = sample } { print }' \\
     read_fraction.tsv > ${meta.id}.read_fraction.tsv
    """
}
