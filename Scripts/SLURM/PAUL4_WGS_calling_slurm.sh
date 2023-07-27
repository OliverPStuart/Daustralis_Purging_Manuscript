#!/bin/bash
#SBATCH --job-name=PAUL4_calling
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

### Calling variants with freebayes on PAUL4_WGS alignment
### We could paralellise but I'm too lazy, and it's not necessary

HOME_DIR=/mnt/data/dayhoff/home/u6905905
WORKING_DIR=${HOME_DIR}/Analyses/PAUL4_WGS
ALN_DIR=${HOME_DIR}/Alignments/WGS
REF_DIR=${HOME_DIR}/References

if [ ! -d "${WORKING_DIR}" ]
then
  mkdir ${WORKING_DIR}
fi

cd ${WORKING_DIR}

${HOME_DIR}/freebayes \
-f ${REF_DIR}/LHISI_Scaffold_Assembly.fasta \
--use-best-n-alleles 3 \
-C 5 \
-0 \
--report-monomorphic \
${ALN_DIR}/PAUL4_WGS.bam > PAUL4_Variants_Raw.vcf
