#!/bin/bash
#SBATCH --job-name=repeat_modeler_LHISI
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=256G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

### Source .bashrc for conda environments
#source ${HOME_DIR}/.bashrc
eval "$(conda shell.bash hook)"

# Find and run script
cd /mnt/data/dayhoff/home/u6905905/Daus_WGS_Paper/Scripts
sh repeat_mask_pipeline.sh
