process BINETTE {
    // Bin refinement (as in the CoMS workflow): combines the bin sets from
    // all binners into one non-redundant set per sample, scoring candidate
    // bins with CheckM2. Final bins are named <sample>_binN.fa, so they can
    // be pooled across samples for dRep without renaming.
    tag "$meta.id"
    label 'process_high'
    publishDir { "${params.outdir}/02_binette/${meta.id}" }, mode: 'copy'

    container "/fs/ess/community/CoMS_microbiome_tools/sif-apps_staging/binette/1.2.1/binette_v1.2.1.sif"

    input:
    tuple val(meta), path(assembly), path(bin_dirs)
    path checkm2_db

    output:
    tuple val(meta), path("${meta.id}_final_bins"), emit: bins
    tuple val(meta), path("${meta.id}_final_bins_quality.tsv"), emit: quality
    tuple val(meta), path("${meta.id}_input_bins_quality"), emit: input_quality

    script:
    """
    # A binner can finish with zero bins; only bin sets that contain bins
    # are given to Binette.
    BIN_DIRS=""
    for d in ${bin_dirs}; do
        if ls "\$d"/*.fa > /dev/null 2>&1; then
            BIN_DIRS="\$BIN_DIRS \$d"
        fi
    done

    if [ -n "\$BIN_DIRS" ]; then
        binette \\
            --bin_dirs \$BIN_DIRS \\
            --contigs ${assembly} \\
            --checkm2_db ${checkm2_db} \\
            --prefix ${meta.id} \\
            --threads ${task.cpus} \\
            --outdir binette_out

        mv binette_out/final_bins ${meta.id}_final_bins
        mv binette_out/final_bins_quality_reports.tsv ${meta.id}_final_bins_quality.tsv
        mv binette_out/input_bins_quality_reports ${meta.id}_input_bins_quality
    else
        echo "No binner produced bins for ${meta.id}; nothing to refine" >&2
        mkdir -p ${meta.id}_final_bins ${meta.id}_input_bins_quality
        printf "name\\tcompleteness\\tcontamination\\n" > ${meta.id}_final_bins_quality.tsv
    fi
    """
}


process MERGE_INPUT_BIN_QUALITY {
    // One table per sample with the CheckM2 quality of every bin from each
    // binner, taken from Binette (which already ran CheckM2 on its inputs).
    tag "$meta.id"
    label 'process_low'
    publishDir { "${params.outdir}/02_binette/${meta.id}" }, mode: 'copy'

    input:
    tuple val(meta), path(reports_dir)

    output:
    tuple val(meta), path("${meta.id}_input_bins_quality.tsv"), emit: table

    script:
    """
    awk -v OFS='\\t' -v sample=${meta.id} '
        FNR == 1 {
            binner = FILENAME
            sub(/.*input_bins_[0-9]+\\./, "", binner); sub(/_bins\\.tsv\$/, "", binner)
            if (NR == 1) print "sample", "binner", \$0
            next
        }
        { print sample, binner, \$0 }
    ' ${reports_dir}/*.tsv > ${meta.id}_input_bins_quality.tsv
    """
}