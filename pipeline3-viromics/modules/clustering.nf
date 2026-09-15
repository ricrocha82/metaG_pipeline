// clustering.nf

// Rapid genome clustering based on pairwise ANI (https://bitbucket.org/berkeleylab/checkv/src/master/)

process CLUSTERING {

  tag "clustering $sample_id"

  conda "bioconda::checkv=1.0.3"
  container "/fs/project/PAS1117/modules/singularity/CheckV-0.8.1-ClusterONLY.sif"

  publishDir params.cluster, mode:'copy', pattern: "*.self-blastn.ani.tsv"
  publishDir params.cluster, mode:'copy', pattern: "*.self-blastn.clusters.fna"
  publishDir params.cluster, mode:'copy', pattern: "*.self-blastn.clusters.tsv"
  publishDir params.cluster, mode:'copy', pattern: "*.self-blastn.tsv"


  input:
    tuple val(sample_id), path(filtered_checkv_fasta)

  output:
    tuple val(sample_id), path ("*.self-blastn.ani.tsv")           , emit: ani_tsv
    tuple val(sample_id), path ("*.self-blastn.clusters.fna")      , emit: cluster_fasta
    tuple val(sample_id), path ("*.self-blastn.clusters.tsv")      , emit: cluster_tsv
    tuple val(sample_id), path ("*.self-blastn.tsv")               , emit: blastn_tsv
    // tuple val(sample_id), path("${sample_id}/combined.fna.gz")   , emit: checkv_fasta

  script:
  """
  CheckV-Deduplication.py -i ${filtered_checkv_fasta} \\
                                  --min-ani 95 \\
                                  --min-tcov 85 \\
                                  --min-qcov 0 \\
                                  -o . \\
                                  -c -f \\
                                  -t $params.threads

  """
}


process CONCAT_FASTA {
  // pool sequences
  tag "combine_all"

  publishDir params.combined, mode: 'copy', pattern: "concatenated_all.fasta"

  input:
    // a list of paths (all the per-sample FASTAs)
    path fasta_list

  output:
    path "concatenated_all.fasta", emit: concat_fasta

  script:
  """
  # Concatenate in a deterministic order
  cat ${fasta_list.join(' ')} > concatenated_all.fasta
  """
}



process MAPPING_LONG_READS {
  tag "MAPPING $sample_id"

  // Pick ONE; container wins if both are set.
  conda "bioconda::minimap=0.2_r124"
  // container '/fs/project/PAS1117/modules/singularity/SAMtools-1.22.sif'  

  input:
    tuple val(sample_id), path(cluster_fasta), path(ref_fasta) 

  output:
    tuple val(sample_id), path("${sample_id}_output.sam"), emit: sam_file

  script:
  """
  /fs/project/PAS1117/bioinformatic_tools/minimap2-2.24_x64-linux/minimap2 -ax map-ont \\
                                      -t $params.threads \\
                                      ${cluster_fasta} \\
                                      ${ref_fasta} > "${sample_id}_output.sam"
  """
}


process MAPPING_SHORT_READS {
  tag "MAPPING $sample_id"

  // Pick ONE; container wins if both are set.
  conda "bioconda::minimap=0.2_r124"
  // container '/fs/project/PAS1117/modules/singularity/SAMtools-1.22.sif' 

  input:
    tuple val(sample_id), path(cluster_fasta), path(ref_fasta1), path(ref_fasta2)

  output:
    tuple val(sample_id), path("${sample_id}_clusters.mmi")
    tuple val(sample_id), path("${sample_id}_output.sam"), emit: sam_file

  script:
  """
    # indexing
    /fs/project/PAS1117/bioinformatic_tools/minimap2-2.24_x64-linux/minimap2 -d "${sample_id}_clusters.mmi" ${cluster_fasta}

    # alignment
    /fs/project/PAS1117/bioinformatic_tools/minimap2-2.24_x64-linux/minimap2 -t $params.threads \\
                                          -N 5 -ax sr "${sample_id}_clusters.mmi" \\
                                          ${ref_fasta1} ${ref_fasta2}> "${sample_id}_output.sam"
  """
}


process SAMTOOLS {

  tag "samtool on $sample_id"

  conda "bioconda::samtools=1.22.1"
  container '/fs/project/PAS1117/modules/singularity/SAMtools-1.22.sif'

  publishDir params.mapping, mode: 'copy', pattern: "${sample_id}_sorted.bam"

  input:
    tuple val(sample_id), path(sam_file)

  output:
    tuple val(sample_id), path("${sample_id}_sorted.bam"), emit: bam_file

  script:
  """
    samtools sort -o "${sample_id}_sorted.bam" "${sample_id}_output.sam" -@ $params.threads
    samtools index "${sample_id}_sorted.bam" -@ $params.threads
  """
}


process COVERM_ABUNDANCE {
    tag "CoverM abundance for $sample_id"

    conda "bioconda::coverm=0.6.1"
    // container 'quay.io/biocontainers/coverm:0.6.1--h07ea13f_2'
    // container "/fs/project/PAS1117/modules/singularity/CoverM-0.6.1.sif"

    publishDir params.abundance, mode: 'copy', pattern: "${sample_id}_coverm_output.tsv"

    input:
        tuple val(sample_id), path(bam_file)

    output:
        tuple val(sample_id), path("${sample_id}_coverm_output.tsv"), emit: coverm_output

    script:
    """
    # Calculate abundance using multiple metrics
    /fs/project/PAS1117/modules/singularity/CoverM-0.6.1.sif contig \\
        -b ${bam_file} \\
        --min-read-percent-identity 0.95 \\
        --min-read-aligned-percent 0.75 \\
        -m mean trimmed_mean covered_bases covered_fraction variance length count reads_per_base rpkm tpm \\
        --output-format sparse \\
        --threads ${params.threads} > "${sample_id}_coverm_output.tsv"
    """
}

// create votu table
process CREATE_ABUND_TB {
  tag "creating abundance table"

  // Pick ONE; container wins if both are set.
  conda "conda-forge::python=3.11 pandas"
  // container 'docker.io/python:3.11-slim'   // (optional if you prefer containers)
  
  publishDir "${params.abundance}/tables", mode: 'copy', pattern: "*.tsv"

  input:
    // val trigger
    path coverm_tables

  output:
    path "Unfiltered_vOTU_table.tsv", emit: unfiltered
    path "Filtered_CF_*_vOTU_table.tsv", emit: filtered

  script:
  """
  mkdir -p coverm_input
    
  # Copy all BAM files to the input directory
  cp ${coverm_tables} coverm_input/

  create_vOTU_table.py \\
                    -i coverm_input \\
                    -o . \\
                    --normalization Trimmed Mean \\
                    --covered_fraction 0.5 
  """
}

  // """
  // create_vOTU_table.py \\
  //                   -i $params.abundance \\
  //                   -o .
  // """