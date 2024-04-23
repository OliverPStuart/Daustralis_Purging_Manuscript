#!/bin/bash
#SBATCH --job-name=repeat_modeler_LHISI
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH --mem=256G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

cd /mnt/data/dayhoff/home/u6905905/Daus_WGS_Paper/
echo "Hello world!" > slurm_test_file
