#!/bin/bash
#SBATCH --job-name=PAUL4_remapping
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

### Mapping PAUL4 WGS resequencing reads
### Single sample, no need to make a Snakefile for this

echo "Setting environment"

### Environment
HOME_DIR=/home/scratch/groups/mikheyev/LHISI
REF_DIR=${HOME_DIR}/References
READ_DIR=${HOME_DIR}/Reads/WGS
ALN_DIR=${HOME_DIR}/Alignments/WGS
CONDA_DIR=/home/u6905905/.conda/envs

echo "Allowing conda environment activation"

### Source .bashrc for conda environments
#source ${HOME_DIR}/.bashrc
eval "$(conda shell.bash hook)"

echo "Making output directory for alignment"

### Make output directory if necessary
if [ ! -d "${ALN_DIR}" ]
then
    mkdir ${ALN_DIR}
fi

echo "Activating conda environment for alignment"

### Activate conda environment
conda activate $CONDA_DIR/aln_env

echo "Making index"

### Make bwa index if necessary
if [ ! -f "${REF_DIR}/LHISI_Scaffold_Assembly_Polished.fasta.bwt" ]
then
  bwa index ${REF_DIR}/LHISI_Scaffold_Assembly_Polished.fasta
fi

echo "Running alignment"

### Align, mark duplicates, sort, and output as bam
bwa mem -t 44 -R "@RG\\tID:PAUL4\\tSM:PAUL4\\tLB:NEXTFLEX\\tPL:ILLUMINA" \
${REF_DIR}/LHISI_Scaffold_Assembly_Polished.fasta ${READ_DIR}/PAUL4_R1.fastq.gz ${READ_DIR}/PAUL4_R2.fastq.gz |
samtools fixmate -m -r -u - - | \
samtools sort -u -@ 44 -T ${ALN_DIR}/PAUL4_WGS_1 - | \
samtools markdup -@ 44 -T ${ALN_DIR}/PAUL4_WGS_2 --reference ${REF_DIR}/LHISI_Scaffold_Assembly_Polished.fasta -r -f ${ALN_DIR}/PAUL4_WGS_dupstats.txt - - | \
samtools view -bh > ${ALN_DIR}/PAUL4_WGS.bam

echo "Indexing alignment"

### Index
samtools index -@ 44 ${ALN_DIR}/PAUL4_WGS.bam

conda deactivate

echo "Activating mosdepth env"

### Now analyse depth
if [ ! -d "${HOME_DIR}/Analysis/AnalysingWGSCoverage" ]
then
  mkdir  ${HOME_DIR}/Analysis/AnalysingWGSCoverage
fi

cd ${HOME_DIR}/Analysis/AnalysingWGSCoverage

conda activate mosdepth_env

echo "Analysing depth"

mosdepth -t 44 -n -Q 10 -b 10000 PAUL4_WGS ${ALN_DIR}/PAUL4_WGS.bam

conda deactivate
