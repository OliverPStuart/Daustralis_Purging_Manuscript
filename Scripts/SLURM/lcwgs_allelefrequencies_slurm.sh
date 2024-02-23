#!/bin/bash
#SBATCH --job-name=lcwgs_allele_frequencies
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=42
#SBATCH --time=72:00:00
#SBATCH --mem=96G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err


### Environment
HOME_DIR=/mnt/data/dayhoff/home/u6905905

### Source .bashrc for conda environments
#source ${HOME_DIR}/.bashrc
eval "$(conda shell.bash hook)"

cd ${HOME_DIR}

conda activate $HOME_DIR/.conda/envs/paul4_env

snakemake -s $HOME_DIR/Daus_WGS_Paper/Scripts/Snakefiles/SnakefileAlleleFrequencies --cores 48 --unlock
snakemake -s $HOME_DIR/Daus_WGS_Paper/Scripts/Snakefiles/SnakefileAlleleFrequencies --cores 48
