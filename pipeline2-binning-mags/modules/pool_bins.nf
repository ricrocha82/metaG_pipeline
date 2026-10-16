process POOL_BINS {
    // Collects every sample's Binette bins into one flat directory for dRep.
    // Bin names already carry the sample prefix (Binette --prefix), so there
    // is no collision risk pooling across samples.
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

process POOL_GENOMEINFO {
    // Builds the genomeInfo.csv dRep needs (genome,completeness,contamination)
    // from Binette's quality reports, which are CheckM2 scores. Columns are
    // looked up by header name; "genome" must match pooled_bins/*.fa exactly.
    label 'process_low'
    publishDir "${params.outdir}/03_drep", mode: 'copy'

    input:
    path quality_reports

    output:
    path "genomeInfo.csv", emit: genome_info

    script:
    """
    echo "genome,completeness,contamination" > genomeInfo.csv
    awk -F'\\t' -v OFS=',' '
        FNR == 1 { for (i = 1; i <= NF; i++) col[\$i] = i; next }
        { print \$col["name"] ".fa", \$col["completeness"], \$col["contamination"] }
    ' ${quality_reports} >> genomeInfo.csv
    """
}
