process HOST_PHIX_SUMMARY {
    // One cross-sample table from the Bowtie2 logs: sample, pct_host, pct_phix.
    // Values are Bowtie2's "overall alignment rate". pct_phix is relative to
    // the host-removed reads, not to the raw reads.
    label 'process_low'
    publishDir "${params.outdir}/00_summary", mode: 'copy'

    input:
    path host_logs
    path phix_logs

    output:
    path "host_phix_summary.tsv", emit: summary

    script:
    """
    echo -e "sample\\tpct_host_reads\\tpct_phix_reads" > host_phix_summary.tsv

    for host_log in *_host.bowtie2.log; do
        sample=\$(basename "\$host_log" _host.bowtie2.log)
        pct_host=\$(grep "overall alignment rate" "\$host_log" | cut -d'%' -f1)
        pct_phix=\$(grep "overall alignment rate" "\${sample}_phix.bowtie2.log" | cut -d'%' -f1)
        echo -e "\${sample}\\t\${pct_host}\\t\${pct_phix}" >> host_phix_summary.tsv
    done
    """
}
