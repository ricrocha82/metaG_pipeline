// process METAWRAP_REFINE {
//     // Bin refinement only -- consolidates the 3 bin sets (MetaBAT2,
//     // COMEBin, SemiBin2) into one refined, non-redundant set per sample
//     // using CheckM completeness/contamination scoring. Takes bin FASTA
//     // directories directly, no conversion step needed (unlike DAS_Tool,
//     // which needed contig2bin TSVs).
//     //
//     // No local .sif for MetaWRAP exists in the PAS1117 module library --
//     // using conda. Real caveat: MetaWRAP has a heavy, somewhat dated
//     // conda dependency tree (its own bundled CheckM 1.0.x database
//     // expectations) -- budget real time for the environment to build
//     // cleanly on first run.
//     //
//     // MetaWRAP's bin_refinement is capped at exactly 3 bin sets (-A/-B/-C)
//     // -- fine here since there are only 3 binners (MetaBAT2/COMEBin/SemiBin2).
//     //
//     // Output directory is named uniquely per sample
//     // (${meta.id}_refined_bins) -- NOT a generic fixed name -- since
//     // pooling across samples later (POOL_BINS) would otherwise hit a
//     // file-staging collision if every sample's output shared the same
//     // literal directory name (this exact bug happened with DAS_Tool
//     // before the fix).
//     tag "$meta.id"
//     label 'process_high'
//     publishDir "${params.binning_outdir}/02_bin_refinement/${meta.id}", mode: 'copy'

//     container "/fs/project/PAS1117/modules/singularity/metawrap-1.3.2-checkm2-fork.sif"

//     containerOptions "--bind ${params.checkm2_db}:/checkm2_db:ro"

//     input:
//     tuple val(meta), path(metabat2_dir), path(comebin_dir), path(semibin2_dir)

//     output:
//     tuple val(meta), path("${meta.id}_refined_bins"), emit: bins
//     tuple val(meta), path("${meta.id}_refine_summary.tsv"), emit: summary

//     script:
//     """

//     BIN_ARGS="-A ${metabat2_dir}"
//     if [ -n "\$(ls -A ${comebin_dir} 2>/dev/null)" ]; then
//         BIN_ARGS="\$BIN_ARGS -B ${comebin_dir}"
//     fi
//     BIN_ARGS="\$BIN_ARGS -C ${semibin2_dir}"

//     metawrap bin_refinement \\
//         -o refine_out \\
//         -t ${task.cpus} \\
//         \$BIN_ARGS \\
//         -c 50 -x 10 \\
//         --keep-ambiguous || true

//     mkdir -p ${meta.id}_refined_bins

//     if [ -d refine_out/metawrap_50_10_bins ] && [ -n "\$(ls -A refine_out/metawrap_50_10_bins/*.fa 2>/dev/null)" ]; then
//         for f in refine_out/metawrap_50_10_bins/*.fa; do
//             cp "\$f" "${meta.id}_refined_bins/${meta.id}_\$(basename \$f)"
//         done
//         cp refine_out/metawrap_50_10_bins.stats ${meta.id}_refine_summary.tsv
//     else
//         echo -e "bin\\tcompleteness\\tcontamination\\tGC\\tlineage\\tN50\\tsize\\tbinner" > ${meta.id}_refine_summary.tsv
//         echo "No bins passed -c 50 -x 10 threshold for ${meta.id}" >&2
//     fi
//     """
// }


process METAWRAP_REFINE {
    // Bin refinement only -- consolidates the 3 bin sets (MetaBAT2,
    // COMEBin, SemiBin2) into one refined, non-redundant set per sample
    // using CheckM completeness/contamination scoring. Takes bin FASTA
    // directories directly, no conversion step needed (unlike DAS_Tool,
    // which needed contig2bin TSVs).
    //
    // No local .sif for MetaWRAP exists in the PAS1117 module library --
    // using conda. Real caveat: MetaWRAP has a heavy, somewhat dated
    // conda dependency tree (its own bundled CheckM 1.0.x database
    // expectations) -- budget real time for the environment to build
    // cleanly on first run.
    //
    // MetaWRAP's bin_refinement is capped at exactly 3 bin sets (-A/-B/-C)
    // -- fine here since there are only 3 binners (MetaBAT2/COMEBin/SemiBin2).
    //
    // Output directory is named uniquely per sample
    // (${meta.id}_refined_bins) -- NOT a generic fixed name -- since
    // pooling across samples later (POOL_BINS) would otherwise hit a
    // file-staging collision if every sample's output shared the same
    // literal directory name (this exact bug happened with DAS_Tool
    // before the fix).
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.binning_outdir}/02_bin_refinement/${meta.id}", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/metawrap-1.3.2-checkm2-fork.sif"

    containerOptions "--bind ${params.checkm2_db}:/checkm2_db:ro"

    input:
    tuple val(meta), path(metabat2_dir), path(comebin_dir), path(semibin2_dir)

    output:
    tuple val(meta), path("${meta.id}_refined_bins"), emit: bins
    tuple val(meta), path("${meta.id}_refine_summary.tsv"), emit: summary

    script:
    """

    BIN_ARGS="-A ${metabat2_dir}"
    NEXT_LETTER="B"
    if [ -n "\$(ls -A ${comebin_dir} 2>/dev/null)" ]; then
        BIN_ARGS="\$BIN_ARGS -B ${comebin_dir}"
        NEXT_LETTER="C"
    fi
    BIN_ARGS="\$BIN_ARGS -\${NEXT_LETTER} ${semibin2_dir}"

    metawrap bin_refinement \\
        -o refine_out \\
        -t ${task.cpus} \\
        \$BIN_ARGS \\
        -c 50 -x 10 \\
        --keep-ambiguous || true

    mkdir -p ${meta.id}_refined_bins

    if [ -d refine_out/metawrap_50_10_bins ] && [ -n "\$(ls -A refine_out/metawrap_50_10_bins/*.fa 2>/dev/null)" ]; then
        for f in refine_out/metawrap_50_10_bins/*.fa; do
            cp "\$f" "${meta.id}_refined_bins/${meta.id}_\$(basename \$f)"
        done
        cp refine_out/metawrap_50_10_bins.stats ${meta.id}_refine_summary.tsv
    else
        echo -e "bin\\tcompleteness\\tcontamination\\tGC\\tlineage\\tN50\\tsize\\tbinner" > ${meta.id}_refine_summary.tsv
        echo "No bins passed -c 50 -x 10 threshold for ${meta.id}" >&2
    fi
    """
}