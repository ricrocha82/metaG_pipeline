#!/bin/bash

#SBATCH --account=PAS1117
#SBATCH --job-name=build_image
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=/fs/project/PAS1117/ricardo/pipeline2-binning-mags/assets/%x_%j.err
#SBATCH --output=/fs/project/PAS1117/ricardo/pipeline2-binning-mags/assets/%x_%j.out




apptainer build /fs/project/PAS1117/modules/singularity/metawrap-1.3.2-checkm2-fork.sif /fs/project/PAS1117/ricardo/pipeline2-binning-mags/assets/metawrap-1.3.2-checkm2-fork.def