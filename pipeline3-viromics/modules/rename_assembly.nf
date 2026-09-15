// rename_assembly.nf
// Pig-wound-specific module — does not touch virus_detection.nf.

process RENAME_ASSEMBLY {

  tag "$sample_id"

  // Reuse the seqkit container already available on OSC — no new image pull needed;
  // `ln` is just standard coreutils, present in any base image.
  container '/fs/project/PAS1117/modules/singularity/seqkit-2.10.0.sif'

  input:
    tuple val(sample_id), path(assembly)   // (sample_id, assembly fasta with arbitrary upstream naming)

  output:
    tuple val(sample_id), path("${sample_id}.fasta"), emit: renamed_fasta   // (sample_id, assembly renamed to match sample_id)

  script:
  """
    # Symlink the assembly to a sample_id-based filename so downstream tools
    # that derive output names from the input file's stem (e.g. geNomad)
    # always produce sample_id-consistent paths, regardless of the
    # upstream assembler's naming convention (e.g. <sample_id>.contigs.fa)
    ln -s ${assembly} ${sample_id}.fasta
  """
}
