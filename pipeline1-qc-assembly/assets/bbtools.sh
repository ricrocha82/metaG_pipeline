#!/bin/bash
#SBATCH --job-name=test_host_clean_reads_v2
#SBATCH --account=PAS1117
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --error=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/code/metaG/%x_%j.err
#SBATCH --output=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/code/metaG/%x_%j.out

# apptainer exec /fs/project/PAS1117/modules/singularity/BBTools-39.31.sif bbmap.sh threads=40 \
#         minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
#         overwrite=t \
#         ref=/fs/project/PAS1117/Pig_burn_july_2025/MSullivan_009_alfonso/host_pig/GCF_000003025.6_Sscrofa11.1_genomic.fna \
#         in1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/assembly/04a_clean_reads/b_4S_clean_1.fastq.gz \
#         in2=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/assembly/04a_clean_reads/b_4S_clean_2.fastq.gz \
#         outu1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_removed_1.fastq.gz \
#         outu1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_removed_2.fastq.gz \
#         outm1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_matched_1.fastq.gz \
#         outm2=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_matched_2.fastq.gz \
#         statsfile=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/hostRemoval.stats \
#         nodisk


apptainer exec /fs/project/PAS1117/modules/singularity/BBTools-39.31.sif bbmap.sh threads=40 \
        minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
        overwrite=t \
        ref=/fs/project/PAS1117/Pig_burn_july_2025/MSullivan_009_alfonso/host_pig/GCF_000003025.6_Sscrofa11.1_genomic.fna \
        in1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_removed_1.fastq.gz \
        in2=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_removed_2.fastq.gz \
        outu1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_removed_1_version2.fastq.gz \
        outu1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_removed_2_version2.fastq.gz \
        outm1=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_matched_1_version2.fastq.gz \
        outm2=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/host_matched_2_version2.fastq.gz \
        statsfile=/fs/project/PAS1117/ricardo/pig_burn/host_depletion/results/hostRemoval_version2.stats \
        nodisk
