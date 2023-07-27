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

### Environment
HOME_DIR=/mnt/data/dayhoff/home/u6905905
REF_DIR=${HOME_DIR}/References
READ_DIR=${HOME_DIR}/Reads/WGS
ALN_DIR=${HOME_DIR}/Alignments/WGS

### Source .bashrc for conda environments
#source ${HOME_DIR}/.bashrc
eval "$(conda shell.bash hook)"

### Make output directory if necessary
if [ ! -d "${ALN_DIR}" ]
then
    mkdir ${ALN_DIR}
fi

### Activate conda environment
conda activate $HOME_DIR/.conda/envs/aln_env

### Make bwa index if necessary
if [ ! -f "${REF_DIR}/LHISI_Scaffold_Assembly.fasta.bwt" ]
then
  bwa index ${REF_DIR}/LHISI_Scaffold_Assembly.fasta
fi

### Align, mark duplicates, sort, and output as bam
bwa mem -t 44 -R "@RG\\tID:PAUL4\\tSM:PAUL4\\tLB:NEXTFLEX\\tPL:ILLUMINA" \
${REF_DIR}/LHISI_Scaffold_Assembly.fasta ${READ_DIR}/PAUL4_R1.fastq.gz ${READ_DIR}/PAUL4_R2.fastq.gz |
samtools fixmate -m -r -u - - | \
samtools sort -u -@ 44 -T ${ALN_DIR}/PAUL4_1 - | \
samtools markdup -@ 44 -T ${ALN_DIR}/PAUL4_2 --reference ${REF_DIR}/LHISI_Scaffold_Assembly.fasta -r -f ${ALN_DIR}/PAUL4_dupstats.txt - - | \
samtools view -bh > ${ALN_DIR}/PAUL4.bam

### Rename PAUL4 high coverage alignment after read tags added
### We want the same read tags as lc data, but different file name

mv PAUL4.bam PAUL4_WGS.bam

### Index
samtools index -@ 44 ${ALN_DIR}/PAUL4.bam

conda deactivate

### Now analyse depth
if [ ! -d "${HOME_DIR}/Analyses/AnalysingWGSCoverage" ]
then
  mkdir  ${HOME_DIR}/Analyses/AnalysingWGSCoverage
fi

cd ${HOME_DIR}/Analyses/AnalysingWGSCoverage

conda activate mosdepth_env

mosdepth -t 44 -n -Q 10 -b 10000 PAUL4_WGS ${ALN_DIR}/PAUL4_WGS.bam
