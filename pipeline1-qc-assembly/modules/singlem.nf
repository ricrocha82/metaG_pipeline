process SINGLEM_PIPE {
    // SingleM-0.13.2.sif is an old version -- confirmed via --help that it
    // has no prokaryotic_fraction/microbial_fraction subcommand at all
    // (not a flag difference, the feature doesn't exist in this version),
    // and uses --otu_table (underscore) not --otu-table. This is now the
    // full SingleM step for this pipeline: an OTU table only, no
    // fraction estimate.
    //
    // Still unverified: whether this container has a default metapackage
    // configured (no --metapackage flag passed here). If this fails with
    // a database/metapackage-not-found error, that's the next thing to
    // diagnose -- see .command.err for the real message.
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}/05b_singlem/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/SingleM-0.13.2.sif"

    input:
    tuple val(meta), path(reads1), path(reads2)

    output:
    tuple val(meta), path("${meta.id}.otu_table.tsv"), emit: otu_table

    script:
    """
    singlem pipe \\
        --forward ${reads1} \\
        --reverse ${reads2} \\
        --otu_table ${meta.id}.otu_table.tsv \\
        --threads ${task.cpus}
    """
}




// process SINGLEM_PIPE {
//     // No --metapackage passed -- relies on SingleM's own default set
//     // (whatever SingleM-0.13.2.sif resolves via SINGLEM_METAPACKAGE_PATH
//     // or its own bundled default). If this errors with something like
//     // "no metapackage found", that means this container does NOT have a
//     // default configured, and a --metapackage path will need adding back.
//     tag "$meta.id"
//     label 'process_high'
//     publishDir "${params.outdir}/05b_singlem/${meta.id}", mode: 'copy'

//     container "/fs/project/PAS1117/modules/singularity/SingleM-0.13.2.sif"

//     input:
//     tuple val(meta), path(reads1), path(reads2)

//     output:
//     tuple val(meta), path(reads1), path(reads2), path("${meta.id}.profile"), emit: profile

//     script:
//     """
//     singlem pipe \\
//         --forward ${reads1} \\
//         --reverse ${reads2} \\
//         --threads ${task.cpus} \\
//         -p ${meta.id}.profile
//     """
// }

// process SINGLEM_PROKARYOTIC_FRACTION {
//     // https://wwood.github.io/singlem/tools/prokaryotic_fraction
//     tag "$meta.id"
//     label 'process_medium'
//     publishDir "${params.outdir}/05b_singlem/${meta.id}", mode: 'copy'

//     container "/fs/project/PAS1117/modules/singularity/SingleM-0.13.2.sif"

//     input:
//     tuple val(meta), path(reads1), path(reads2), path(profile)

//     output:
//     tuple val(meta), path("${meta.id}.spf.tsv"), emit: fraction

//     script:
//     """
//     singlem prokaryotic_fraction \\
//         --forward ${reads1} \\
//         --reverse ${reads2} \\
//         -p ${profile} \\
//         --output-tsv ${meta.id}.spf.tsv
//     """
// }
