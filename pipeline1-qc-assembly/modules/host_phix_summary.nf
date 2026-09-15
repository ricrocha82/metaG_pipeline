process HOST_PHIX_SUMMARY {
    // Parses the BBMap/BBDuk stats files into one clean cross-sample table:
    // sample, pct_host, pct_phix
    //
    // CAVEAT: I don't have BBTools available to test-run and confirm exact
    // stats file column layouts, so this parsing is best-effort based on
    // BBMap/BBDuk's documented output format, not verified against real
    // output. After your first real run, open one *_hostRemoval.stats and
    // one *_PhiXFiltering.stats file and sanity-check the numbers in
    // host_phix_summary.tsv actually match what's in those raw files --
    // if the awk/grep parsing is off, the raw per-sample .stats files
    // (always published regardless) are the ground truth to fall back on.
    label 'process_low'
    publishDir "${params.outdir}/00_summary", mode: 'copy'

    input:
    path host_stats_files
    path phix_stats_files

    output:
    path "host_phix_summary.tsv", emit: summary

    script:
    """
    echo -e "sample\\tpct_host_reads\\tpct_phix_reads" > host_phix_summary.tsv

    for f in ${host_stats_files}; do
        sample=\$(basename \$f _hostRemoval.stats)
        # "mapped:" line under "Read 1 data" = % of reads that aligned to
        # the host reference, i.e. host contamination %
        pct_host=\$(grep -A3 "^Read 1 data" \$f | grep "^mapped:" | awk '{print \$2}' | tr -d '%')
        pct_host=\${pct_host:-NA}

        phix_file="\$(echo \${f%_hostRemoval.stats})_PhiXFiltering.stats"
        # match against the corresponding PhiX stats file for this sample
        for pf in ${phix_stats_files}; do
            pf_sample=\$(basename \$pf _PhiXFiltering.stats)
            if [ "\$pf_sample" == "\$sample" ]; then
                # BBDuk's stats= output: total matched reads on the summary line
                pct_phix=\$(awk -F'\\t' '/^#/{next} {sum+=\$5} END{printf "%.4f", sum}' \$pf)
                pct_phix=\${pct_phix:-NA}
            fi
        done

        echo -e "\${sample}\\t\${pct_host}\\t\${pct_phix}" >> host_phix_summary.tsv
    done
    """
}
