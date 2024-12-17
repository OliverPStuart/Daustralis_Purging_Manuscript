#!/bin/bash
#SBATCH --job-name=sitewise_h_all
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Run after alignments and per-individual genomic H estimate is complete

######################################
####          ENV SETUP           ####
######################################

HOME_DIR=/mnt/data/dayhoff/home/scratch/groups/mikheyev/LHISI
ALN_DIR=${HOME_DIR}/Alignments/WGS
REF_DIR=${HOME_DIR}/References
ANALYSIS_DIR=${HOME_DIR}/Analysis
SOFT_DIR=/mnt/data/dayhoff/home/u6905905/
ANGSD=${SOFT_DIR}/angsd
REFERENCE=${REF_DIR}/LHISI_Scaffold_Assembly.fasta

conda init bash
eval "$(conda shell.bash hook)"
conda activate ${SOFT_DIR}/.conda/envs/paul4_env

mkdir -p ${ANALYSIS_DIR}/SitewiseH
cd ${ANALYSIS_DIR}/SitewiseH

######################################
####     PART 1: Per-site H       ####
######################################

# Make a file containing the names of sample bams
# Only include samples we know ran well

find $ALN_DIR | grep bam$ | grep "C01230\|C01225\|C10223\|C01211\|C01217\|C01220\|C01234\|C01224" > lhip_bam_list.txt
find $ALN_DIR | grep bam$ | grep "C01210\|C01215\|C01223\|C01222\|C01232\|C01231\|C01219\|C01227\|C01213" > lhisi_bam_list.txt

# Make bed file of sites
# This is all sites in the genome, but excluding repetitive and non-mappable regions

head -n 16 ${REF_DIR}/LHISI_Scaffold_Assembly.fasta.fai | awk '{print $1"\t"1"\t"$2}' > autosomes.bed
bedtools subtract -a autosomes.bed -b ${REF_DIR}/problematic_regions.bed > tmp
bedtools subtract -a tmp -b ${REF_DIR}/Annotation/repeats_major.bed | \
awk '{print $1 "\t" $2+1 "\t" $3+1}' > good_sites
rm tmp

# Weird access time overwriting in ANGSD, use sleep to avoid error message

${ANGSD}/angsd sites index good_sites
sleep 5s
touch good_sites.*

# Now run the analysis separately for both populations

for pop in lhisi lhip
do
${ANGSD}/angsd -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 30 -minQ 30 \
 -b ${pop}_bam_list.txt -P 12 -out ${pop}_hwe -ref ${REFERENCE} -anc ${REFERENCE} -sites good_sites \
 -doCounts 1 -setMinDepthInd 3 -setMaxDepthInd 8 -doMajorMinor 2 -doHWE 1 -maxHetFreq 1 -GL 2 2> ${pop}_hwe.err
done

######################################
## PART 2: Per-ind H by site class  ##
######################################

# Extract depth files, requires per individual H estimation to be complete
mkdir -p depth_files
ARCHIVE=$(find ${ANALYSIS_DIR}/SingleSampleHEstimation -name 'SingleSampleH_*.tar.gz' -type f | head -n 1)
tar -xzvf ${ARCHIVE} -C depth_files --wildcards '*_nonRepeat'
# Lots of files, so it takes a bit

# Get files for coordinates of all site types
TYPE_ARR=(HIGH MODERATE LOW MODIFIER)
for TYPE in "${TYPE_ARR[@]}"
do
  awk -v TYPE=${TYPE} 'NR >1 && $4 == TYPE {print $1"\t"$2-1"\t"$2}' ${ANALYSIS_DIR}/DeleteriousMutations/vars_classified.txt > ${TYPE}.bed
done

# Make a regions file to restrict angsd analysis
CHROM_ARR=(CM056993.1 CM056994.1 CM056995.1 CM056996.1 CM056997.1 CM056998.1 CM056999.1 CM057000.1 CM057001.1 CM057002.1 CM057003.1 CM057004.1 CM057005.1 CM057006.1 CM057007.1 CM057008.1)
for CHROM in "${CHROM_ARR[@]}"; do
echo "${CHROM}:1-" >> regions_file
done

# Now loop over all individuals, depths, chromosomes, and site types
IND_ARR=(C01210 C01211 C01213 C01215 C01217 C01219 C01220 C01222 C01223 C01224 C01225 C01227 C01230 C01231 C01232 C01234 C10223 PAUL4 VAN2)
DEPTH_ARR=(low mid)

for IND in "${IND_ARR[@]}"; do
  for DEPTH in "${DEPTH_ARR[@]}"; do
    for TYPE in "${TYPE_ARR[@]}"; do

      #Prefix to save time
      PREFIX=${IND}_${DEPTH}_${TYPE}

      # Loop over all chroms and paste sites into bed file
      for CHROM in "${CHROM_ARR[@]}"; do

        # Get overlapping sites at chrom, depth, and type
        bedtools intersect -a ${TYPE}.bed -b depth_files/${IND}_${DEPTH}Depth_${CHROM}_nonRepeat >> ${PREFIX}.bed

      done

      # Index the sites
      ${ANGSD}/angsd sites index ${PREFIX}.bed

      # Get sample allele frequencies per site
      ${ANGSD}/angsd \
      -i ${ALN_DIR}/${IND}.bam \
      -anc ${REFERENCE} -ref ${REFERENCE} \
      -out ${PREFIX} \
      -nThreads 8  \
      -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 30 -minQ 30 \
      -doSaf 1 -GL 2 -rf regions_file -sites ${PREFIX}.bed 2> ${PREFIX}_saf.err

      # Get single sample SFS
      ${ANGSD}/misc/realSFS \
      -P 8 -fold 1 -maxIter 5000 \
      ${PREFIX}.saf.idx > ${PREFIX}.ml 2> ${PREFIX}_ml.err

      # Collate info into table
      FILECONTENTS=$(<"${PREFIX}.ml")
      echo -e "${IND}\t${DEPTH}\t${TYPE}\t${FILECONTENTS}" >> site_hets.txt

    done ; done ; done
