process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/00_multiqc", mode: 'copy'

    container "/fs/project/PAS1117/modules/singularity/MultiQC-1.7.sif"

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc .
    """
}
