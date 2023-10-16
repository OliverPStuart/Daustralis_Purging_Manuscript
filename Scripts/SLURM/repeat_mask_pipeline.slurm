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

### Evironment
cd /mnt/data/dayhoff/home/u6905905/Daus_WGS_Paper
. ../shell.config
cd ${SCRATCH_DIR}/References/Annotation

# Download singularity container

singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest

# Now run repeatmodeler

singularity exec repeatmodeler.sif BuildDatabase -name LHISI_Scaffold_Assembly_Polished
singularity exec repeatmodeler.sif RepeatModeler -database LHISI_Scaffold_Assembly_Polished -threads 28

# Now use CD-hit to reduce redundancy with the previous run

conda activate cdhit_env
cat LHISI_Scaffold_Assembly_Polished-families.fa LHISI_Scaffold_Assembly.repeats_filtered.fasta > combined_repeats.fa
cd-hit-est combined_repeats.fa -c 0.99 -o combined_repeats_clustered.fa
conda deactivate
