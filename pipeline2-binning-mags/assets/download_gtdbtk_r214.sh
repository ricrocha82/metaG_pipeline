#!/bin/bash
#SBATCH --account=PAS1117
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=download_gtdbtk_r207
#SBATCH --error=/fs/project/PAS1117/ricardo/pipeline2-binning-mags/assets/%x_%j.err
#SBATCH --output=/fs/project/PAS1117/ricardo/pipeline2-binning-mags/assets/%x_%j.out

mkdir -p /fs/ess/PAS1117/ricardo/gtdbtk_r214_data
cd /fs/ess/PAS1117/ricardo/gtdbtk_r214_data
wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
tar xvzf gtdbtk_r214_data.tar.gz

# https://ecogenomics.github.io/GTDBTk/installing/index.html