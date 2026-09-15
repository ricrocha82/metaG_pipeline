
// metapop
process METAPOP_MICRODIVERSITY {
  tag "$method"

  container "/fs/project/PAS1117/ricardo/virion3/conf/metapop/metapop.sif"
  publishDir "${params.microdiversity}/${method}", mode: 'copy'

  input:
  tuple val(method), path(bam_files), path(cluster_fasta)

  output:
  path "MetaPop/10.Microdiversity/*.tsv"

  script:
  """
  mkdir -p bam_input reference_dir
  cp ${bam_files} bam_input/
  cp ${cluster_fasta} reference_dir/
  metapop --input_samples bam_input \\
          --reference reference_dir \\
          --threads $params.threads \\
          --no_macro --no_viz \\
          --output .
  """
}