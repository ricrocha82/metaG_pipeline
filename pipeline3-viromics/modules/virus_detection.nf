// virus_detection.nf

process VIRSORTER {

  tag "Virsorter2 on $sample_id"

  conda "bioconda::virsorter=2.2.4"
  container '/fs/project/PAS1117/modules/singularity/VirSorter2-2.2.3.sif'

  publishDir params.vs2, mode:'copy', pattern: "${sample_id}/final-viral-score.tsv", saveAs: { filename -> "${sample_id}_final-viral-score.tsv" }
  publishDir params.vs2, mode:'copy', pattern: "${sample_id}/final-viral-boundary.tsv", saveAs: { filename -> "${sample_id}_final-viral-boundary.tsv" }
  publishDir params.vs2, mode:'copy', pattern: "${sample_id}/final-viral-combined.fa", saveAs: { filename -> "${sample_id}_final-viral-combined.fa" }

  input:
    tuple val(sample_id), path(read_path)

  output:
    tuple val(sample_id), path("${sample_id}/final-viral-score.tsv") , emit: viral_score
    tuple val(sample_id), path("${sample_id}/final-viral-boundary.tsv") , emit: viral_boundary
    tuple val(sample_id), path("${sample_id}/final-viral-combined.fa") , emit: viral_combined

  script:
  """
  virsorter run --keep-original-seq -i ${read_path} \
                                        -w ${sample_id} \
                                        --include-groups dsDNAphage,ssDNA \
                                        --min-length $params.vs2_length \
                                        --min-score 0.5 \
                                        -j $params.threads all
  """
}

process GENOMAD_END_TO_END {

  tag "geNomad on $sample_id"

  conda "bioconda::genomad=1.9.0"
  container "/fs/project/PAS1117/modules/singularity/GeNomad-1.9.0.sif"

  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_virus_summary.tsv", saveAs: { filename -> "${sample_id}_virus_summary.tsv" }
  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_virus.fna", saveAs: { filename -> "${sample_id}_virus.fna" }
  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_virus_genes.tsv", saveAs: { filename -> "${sample_id}_virus_genes.tsv" }
  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_virus_proteins.faa", saveAs: { filename -> "${sample_id}_virus_proteins.faa" }
  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_plasmid.fna", saveAs: { filename -> "${sample_id}_plasmid.fna" }
  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_plasmid_summary.tsv", saveAs: { filename -> "${sample_id}_plasmid_summary.tsv" }
  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_plasmid_genes.tsv", saveAs: { filename -> "${sample_id}_plasmid_genes.tsv" }
  publishDir params.genomad, mode:'copy', pattern: "${sample_id}/${sample_id}_summary/${sample_id}_plasmid_proteins.faa", saveAs: { filename -> "${sample_id}_plasmid_proteins.faa" }
  

  input:
    tuple val(sample_id), path(read_path)
    val  genomad_db

  output:
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_virus_summary.tsv") , emit:  virus_summary
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_virus.fna") , emit:  genomad_virus
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_plasmid.fna") , emit:  genomad_plasmid
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_plasmid_summary.tsv") , emit:  plasmid_summary
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_virus_genes.tsv"), emit: virus_genes
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_plasmid_genes.tsv"), emit: plasmid_genes
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_virus_proteins.faa"), emit: virus_proteins
    tuple val(sample_id), path("${sample_id}/${sample_id}_summary/${sample_id}_plasmid_proteins.faa"), emit: plasmid_proteins

  script:
  """
    genomad end-to-end \\
        --cleanup \\
        --splits $params.threads \\
        ${read_path} \\
        ${sample_id} \\
        $genomad_db
  """
}

process COMBINE_VS2_GENOMAD {
  tag "combine VS2 + geNomad"

  // Pick ONE; container wins if both are set.
  conda "conda-forge::python=3.11 pandas"
  // container 'docker.io/python:3.11-slim'   // (optional if you prefer containers)
  
  publishDir params.combined, mode: 'copy', pattern: "*_VS2_GENOMAD.txt"

  input:
    tuple val(sample_id), path(virus_score), path(virus_summary)

  output:
    path "*_VS2_GENOMAD.txt", emit: combined_df
    // stdout emit: my_stdout_output // Emit stdout

  script:
  """
  echo ${virus_score}
  echo ${virus_summary}
  combine_virus_id.py \\
            --vs2-score ${virus_score} \\
            --genomad-summary ${virus_summary} \\
            --samples $params.input_csv \\
            --outdir .
  """
}

// """
//   combine_virus_id.py \\
//             --vs2-dir $params.vs2 \\
//             --genomad-dir $params.genomad \\
//             --samples $params.input_csv \\
//             --outdir .

//     """


process SELECT_CONTIGS {
  tag "$sample_id"

  // Pick ONE; container wins if both are set.
  conda "conda-forge::python=3.11 pandas"  

  // publishDir params.combined, mode: 'copy', pattern: "${sample_id}_VS2_GENOMAD.fasta"
  publishDir params.combined, mode: 'copy', pattern: "${sample_id}_contigs.txt"

  input:
    tuple val(sample_id), path(combined_df), path(read_path)

  output:
    tuple val(sample_id), path("${sample_id}_contigs.txt"), emit: out_df_filtered
    // tuple val(sample_id), path("${sample_id}_VS2_GENOMAD.fasta"), emit: to_concat

  script:
  """
  filter_vs2_genomad.py --input ${combined_df} --sample ${sample_id}
  """
}


process SELECT_CONTIGS_GENOMAD {
  tag "$sample_id"

  // Pick ONE; container wins if both are set.
  conda "conda-forge::python=3.11 pandas"  

  // publishDir params.combined, mode: 'copy', pattern: "${sample_id}_VS2_GENOMAD.fasta"
  publishDir params.combined, mode: 'copy', pattern: "${sample_id}_contigs.txt"

  input:
    tuple val(sample_id), path(virus_summary), path(read_path)

  output:
    tuple val(sample_id), path("${sample_id}_contigs.txt"), emit: out_df_filtered
    // tuple val(sample_id), path("${sample_id}_VS2_GENOMAD.fasta"), emit: to_concat

  script:
  """
  filter_genomad.py --input ${virus_summary} \\
                    --sample ${sample_id} \\
                    --length $params.virus_length \\
                    --virus-score $params.virus_score
  """
}


process EXTRACT_CONTIGS {
  tag "$sample_id"
  container '/fs/project/PAS1117/modules/singularity/seqkit-2.10.0.sif'

  publishDir params.combined, mode: 'copy', pattern: "${sample_id}_VS2_GENOMAD.fasta"

  input:
    tuple val(sample_id), path(out_df_filtered), path(read_path), path(genomad_virus)

  output:
    tuple val(sample_id), path("${sample_id}_VS2_GENOMAD.fasta"), emit: filt_fasta

  script:
  """
  seqkit grep -f "${out_df_filtered}" "${genomad_virus}" > "${sample_id}_VS2_GENOMAD.fasta"
  """
}

// process FILTER_VS2_GENOMAD {
//   tag "$sample_id"

//   // Pick ONE; container wins if both are set.
//   conda "bioconda::seqkit=2.10.1"
//   container '/fs/project/PAS1117/modules/singularity/seqkit-2.10.0.sif'  

//   publishDir params.combined, mode: 'copy', pattern: "${sample_id}_VS2_GENOMAD.fasta"
//   publishDir params.combined, mode: 'copy', pattern: "${sample_id}_contigs.txt"

//   input:
//     tuple val(sample_id), path(combined_df), path(read_path)

//   output:
//     tuple val(sample_id), path("${sample_id}_contigs.txt"), emit: out_df_filtered
//     tuple val(sample_id), path("${sample_id}_VS2_GENOMAD.fasta"), emit: to_concat

//   script:
//   """
//   # -F'\t' for TSV files
//   awk '{ if ((\$2 >= 0.5) && (\$7 >= 5000) && (\$9 >= 0.7)) { print \$1}}' ${combined_df} > ${sample_id}_contigs.txt

//   seqkit grep -nrif ${sample_id}_contigs.txt ${read_path} -o ${sample_id}_VS2_GENOMAD.fasta
//   """
// }
