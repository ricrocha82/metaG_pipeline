process METABAT2 {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.binning_outdir}/01_binning/${meta.id}/metabat2", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/MetaBAT2-2.14.sif"

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)

    output:
    tuple val(meta), path("metabat2_bins"), emit: bins

    script:
    """
    jgi_summarize_bam_contig_depths --outputDepth depth.txt ${bam}
    mkdir -p metabat2_bins
    metabat2 -i ${assembly} -a depth.txt \\
        -o metabat2_bins/${meta.id}_bin \\
        -t ${task.cpus} \\
        -m ${params.binning_min_contig_size ?: 1500} \\
+       --seed ${params.binning_seed ?: 1}
    """
}

process COUNT_CONTIGS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)

    output:
    tuple val(meta), path(assembly), path(bam), path(bai), path("n_contigs.txt"), emit: counted

    script:
    """
    grep -c '^>' ${assembly} > n_contigs.txt
    """
}


// process COMEBIN {
//     // Conda env, per your request -- no local .sif for COMEBin exists on
//     // PAS1117. Run with -profile conda (see README).
//     errorStrategy 'ignore'
//     tag "$meta.id"
//     label 'process_high'
//     publishDir "${params.binning_outdir}/01_binning/${meta.id}/comebin", mode: 'copy'

//     container "/fs/project/PAS1117/modules/singularity/comebin-1.0.4.sif" 

//     input:
//     tuple val(meta), path(assembly), path(bam), path(bai)

//     output:
//     tuple val(meta), path("comebin_bins"), emit: bins
//     tuple val(meta), path("comebin_quality"), emit: comebin_quality

//     script:
//     """
//     mkdir -p comebin_bins comebin_quality

//     run_comebin.sh -a ${assembly} -o comebin_out -p . -t ${task.cpus} > comebin_run.log 2>&1

//     cp comebin_out/comebin_res/comebin_res_bins/*.fa comebin_bins/ 2>/dev/null || true

//     n_bins=\$(ls comebin_bins/*.fa 2>/dev/null | wc -l)
//     if [ "\$n_bins" -eq 0 ]; then
//         echo "ERROR: COMEBin produced zero bins for ${meta.id} despite passing the contig-count gate -- this is unexpected, check comebin_run.log" >&2
//         find comebin_out -maxdepth 4 -type d >&2
//         exit 1
//     fi
//     echo "COMEBin produced \$n_bins bins"

//     cp comebin_out/comebin_res/cluster_res/unitem_profile/*_quality.tsv comebin_quality/ 2>/dev/null || true
//     cp comebin_out/comebin_res/cluster_res/unitem_profile/bin_quality_summary.tsv comebin_quality/ 2>/dev/null || true
//     cp comebin_out/comebin_res/cluster_res/unitem_profile/binning_methods/*/checkm_bac/genome_quality.tsv \\
//         comebin_quality/checkm_bacteria_genome_quality.tsv 2>/dev/null || true
//     cp comebin_out/comebin_res/cluster_res/unitem_profile/binning_methods/*/checkm_ar/genome_quality.tsv \\
//         comebin_quality/checkm_archaea_genome_quality.tsv 2>/dev/null || true
//     """
// }

process COMEBIN {
    errorStrategy 'ignore'
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.binning_outdir}/01_binning/${meta.id}/comebin", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/comebin-1.0.4.sif"

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)

    output:
    tuple val(meta), path("comebin_bins"), emit: bins
    tuple val(meta), path("comebin_quality"), emit: comebin_quality

    script:
    """
    mkdir -p comebin_bins comebin_quality

    run_comebin.sh -a ${assembly} -o comebin_out -p . -t ${task.cpus} > comebin_run.log 2>&1

    cp comebin_out/comebin_res/comebin_res_bins/*.fa comebin_bins/ 2>/dev/null || true

    n_bins=\$(ls comebin_bins/*.fa 2>/dev/null | wc -l)
    if [ "\$n_bins" -eq 0 ]; then

        # IMPORTANT: never use 'exit 1' here. Nextflow only caches (for
        # -resume) tasks that finish with exit 0 -- a task that fails never
        # becomes a valid cache entry, so with exit 1 COMEBin would be
        # re-run on EVERY future -resume, even with errorStrategy 'ignore'.
        # By treating "zero bins" as a valid result (exit 0), the failure
        # becomes part of the cache normally, and -resume simply reuses
        # the empty result.
        
        if grep -qE "IndexError: list index out of range|max\\(\\) arg is an empty sequence" comebin_run.log; then
            echo "KNOWN BUG: COMEBin marker-gene seeding failure for ${meta.id} (ziyewang/COMEBin#40) -- treating as zero bins, cached, will NOT re-run on resume" >&2
        else
            echo "WARNING: COMEBin produced zero bins for ${meta.id} for an UNRECOGNIZED reason -- check comebin_run.log manually, this is NOT the known marker-gene bug signature" >&2
        fi
        find comebin_out -maxdepth 4 -type d >&2 || true
    else
        echo "COMEBin produced \$n_bins bins"
    fi

    cp comebin_out/comebin_res/cluster_res/unitem_profile/*_quality.tsv comebin_quality/ 2>/dev/null || true
    cp comebin_out/comebin_res/cluster_res/unitem_profile/bin_quality_summary.tsv comebin_quality/ 2>/dev/null || true
    cp comebin_out/comebin_res/cluster_res/unitem_profile/binning_methods/*/checkm_bac/genome_quality.tsv \\
        comebin_quality/checkm_bacteria_genome_quality.tsv 2>/dev/null || true
    cp comebin_out/comebin_res/cluster_res/unitem_profile/binning_methods/*/checkm_ar/genome_quality.tsv \\
        comebin_quality/checkm_archaea_genome_quality.tsv 2>/dev/null || true
    """
}

process SEMIBIN2 {
    // Conda env, per your request -- no local .sif for SemiBin2 exists on
    // PAS1117. Run with -profile conda (see README).
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.binning_outdir}/01_binning/${meta.id}/semibin2", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/semibin2-2.3.0.sif"

    input:
    tuple val(meta), path(assembly), path(bam), path(bai)

    output:
    tuple val(meta), path("semibin2_bins"), emit: bins
    tuple val(meta), path("semibin2_summary"), emit: semibin2_summary

    script:
    """
    SemiBin2 single_easy_bin \\
        -i ${assembly} -b ${bam} \\
        -o semibin2_out \\
        -t ${task.cpus} \\
        --random-seed ${params.binning_seed ?: 1} \\
        --environment global

    mkdir -p semibin2_bins
    cp semibin2_out/output_bins/*.fa.gz semibin2_bins/ 2>/dev/null || true
    gunzip semibin2_bins/*.fa.gz 2>/dev/null || true

    mkdir -p semibin2_summary
    cp semibin2_out/recluster_bins_info.tsv semibin2_summary/ 2>/dev/null || true
    cp semibin2_out/contig_bins.tsv semibin2_summary/ 2>/dev/null || true
    """
}

// process VAMB {
    
//     tag "$meta.id"
//     label 'process_high'
//     publishDir "${params.binning_outdir}/01_binning/${meta.id}/vamb", mode: 'copy'

//     container "/fs/project/PAS1117/modules/singularity/VAMB-3.0.2.sif"

//     input:
//     tuple val(meta), path(assembly), path(bam), path(bai)

//     output:
//     tuple val(meta), path("vamb_bins"), emit: bins

//     script:
//     """
//     mkdir -p bam_dir
//     cp ${bam} ${bai} bam_dir/

//     vamb bin default \\
//         --outdir vamb_out \\
//         --fasta ${assembly} \\
//         -p ${task.cpus}

//     mkdir -p vamb_bins
//     cp vamb_out/bins/*.fna vamb_bins/ 2>/dev/null || true
//     """
// }
