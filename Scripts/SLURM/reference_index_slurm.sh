#!/bin/bash
#SBATCH --job-name=lhisi_reference_index
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

### Environment
HOME_DIR=/mnt/data/dayhoff/home/u6905905

### Source .bashrc for conda environments
#source ${HOME_DIR}/.bashrc
eval "$(conda shell.bash hook)"
conda activate $HOME_DIR/.conda/envs/paul4_env

cd /mnt/data/dayhoff/home/scratch/groups/mikheyev/LHISI/References/

bwa index LHISI_Scaffold_Assembly.fasta
