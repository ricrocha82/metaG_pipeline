process POOL_BINS {
    // Collects every sample's DAS_Tool bins into one flat directory for
    // dRep. Bin filenames already carry their sample prefix (DAS_Tool's
    // -o ${meta.id} flag), so no collision risk pooling across samples.
    label 'process_low'

    input:
    path bin_dirs

    output:
    path "pooled_bins", emit: bins

    script:
    """
    mkdir -p pooled_bins
    for d in ${bin_dirs}; do
        cp \$d/*.fa pooled_bins/ 2>/dev/null || true
    done
    """
}

process POOL_CHECKM2_GENOMEINFO {
    // Reformats all samples' CheckM2 quality_report.tsv files into the
    // single genomeInfo.csv (genome,completeness,contamination) dRep
    // requires, with "genome" values matching pooled_bins/*.fa exactly.
    label 'process_low'
    publishDir "${params.binning_outdir}/03_checkm2", mode: 'copy'

    input:
    path checkm2_results

    output:
    path "genomeInfo.csv", emit: genome_info

    script:
    """
    echo "genome,completeness,contamination" > genomeInfo.csv
    for f in ${checkm2_results}; do
        tail -n +2 "\$f" | awk -F'\\t' -v OFS=',' '{print \$1".fa", \$2, \$3}' >> genomeInfo.csv
    done
    """
}
