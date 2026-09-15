// virus_detection.nf

process CHECKV_END_TO_END {

  tag "checkv $sample_id"

  conda "bioconda::checkv=1.0.3"
  container "/fs/project/PAS1117/modules/singularity/CheckV-1.0.3-updated.sif"

  publishDir params.checkv, mode:'copy', pattern: "quality_summary.tsv", saveAs: { filename -> "${sample_id}_quality_summary.tsv" }
  publishDir params.checkv, mode:'copy', pattern: "completeness.tsv", saveAs: { filename -> "${sample_id}_completeness.tsv" }
  publishDir params.checkv, mode:'copy', pattern: "contamination.tsv", saveAs: { filename -> "${sample_id}_contamination.tsv" }
  publishDir params.checkv, mode:'copy', pattern: "complete_genomes.tsv", saveAs: { filename -> "${sample_id}_complete_genomes.tsv" }
  publishDir params.checkv, mode:'copy', pattern: "combined.fna.gz", saveAs: { filename -> "${sample_id}_combined.fna.gz" }

  input:
    tuple val(sample_id), path(filt_fasta)

  output:
    tuple val(sample_id), path ("quality_summary.tsv")    , emit: quality_summary
    tuple val(sample_id), path ("completeness.tsv")       , emit: completeness
    tuple val(sample_id), path ("contamination.tsv")      , emit: contamination
    tuple val(sample_id), path ("complete_genomes.tsv")   , emit: complete_genomes
    tuple val(sample_id), path ("viruses.fna")            , emit: viruses
    tuple val(sample_id), path ("proviruses.fna")         , emit: proviruses
    tuple val(sample_id), path ("combined.fna.gz")        , emit: checkv_fasta
    // tuple val(sample_id), path("${sample_id}/combined.fna.gz")   , emit: checkv_fasta

  script:
  """
  checkv end_to_end \\
        ${filt_fasta} \\
        . \\
        -t $params.threads

  cat "proviruses.fna" "viruses.fna" > "combined.fna" 

  gzip "combined.fna"
  """
}

process FILTER_CHECKV {
  tag "filtering checkV $sample_id"

  // Pick ONE; container wins if both are set.
  conda "bioconda::seqkit=2.10.1"
  container '/fs/project/PAS1117/modules/singularity/seqkit-2.10.0.sif'  

  publishDir params.checkv, mode: 'copy', pattern: "${sample_id}_checkv_contigs.txt"
  publishDir params.checkv, mode: 'copy', pattern: "${sample_id}_final_viral_contigs.fasta"

  input:
    tuple val(sample_id), path (quality_summary), path (checkv_fasta)

  output:
    tuple val(sample_id), path ("${sample_id}_checkv_contigs.txt") , emit: filtered_checkv_txt
    tuple val(sample_id), path ("${sample_id}_final_viral_contigs.fasta") , emit: filtered_checkv_fasta

  script:
  """
  # -F'\t' for TSV files
  awk -F'\t' '(\$2 >= $params.length) && (\$10 >= $params.virus_completeness) { print \$1 }' ${quality_summary} > "${sample_id}_checkv_contigs.txt"

  seqkit grep -f "${sample_id}_checkv_contigs.txt" ${checkv_fasta} > "${sample_id}_final_viral_contigs.fasta"
  """
}




// process FILTER_CHECKV {
//   tag "filtering checkV"

//   // Pick ONE; container wins if both are set.
//   conda "bioconda::seqkit=2.10.1"
//   container '/fs/project/PAS1117/modules/singularity/seqkit-2.10.0.sif'  

//   publishDir params.checkv, mode: 'copy', pattern: "checkv_contigs.txt"
//   publishDir params.checkv, mode: 'copy', pattern: "final_viral_contigs.fasta"

//   input:
//     path checkv_fasta
//     path quality_summary 

//   output:
//     path "checkv_contigs.txt" , emit: filtered_checkv_txt
//     path "final_viral_contigs.fasta" , emit: filtered_checkv_fasta

//   script:
//   """
//   # -F'\t' for TSV files
//   awk -F'\t' '(\$2 >= 5000) && (\$5 >= 13) && (\$6 >= 1) && (\$7 < 2) && (\$10 >= 70) { print \$1 }' ${quality_summary} > checkv_contigs.txt

//   seqkit grep -f checkv_contigs.txt ${checkv_fasta} > final_viral_contigs.fasta
//   """
// }


