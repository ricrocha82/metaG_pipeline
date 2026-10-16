process SELECT_CANDIDATES {
    // Per sample: decide which sequence represents each detected contig
    // and record which tool found it, keeping every score/annotation column
    // geNomad and VS2 report. Outputs are pooled in main.nf.
    // Python helpers run in the CheckV image (stdlib only).
    tag "$meta.id"
    label 'process_low'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/CheckV/0.8.1/CheckV-0.8.1.sif"

    input:
    tuple val(meta), path(assembly), path(genomad_virus), path(genomad_scores), path(genomad_summary), path(vs2_virus), path(vs2_scores)

    output:
    tuple val(meta), path("${meta.id}_candidates.fna"), emit: fasta
    tuple val(meta), path("${meta.id}_candidates.tsv"), emit: table

    script:
    def genomad_args = genomad_virus ? "--genomad-virus ${genomad_virus} --genomad-scores ${genomad_scores} --genomad-summary ${genomad_summary}" : ''
    def vs2_args     = vs2_virus ? "--vs2-virus ${vs2_virus} --vs2-scores ${vs2_scores}" : ''
    """
    select_candidates.py \\
        --sample ${meta.id} \\
        --contigs ${assembly} \\
        ${genomad_args} \\
        ${vs2_args} \\
        --vs2-groups ${params.vs2_groups} \\
        --out-prefix ${meta.id}_candidates
    """
}
