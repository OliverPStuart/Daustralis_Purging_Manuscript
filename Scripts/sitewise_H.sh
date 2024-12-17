### Estimating site_wise H

# Run after single sample H and deleterious mutation analysis

HOME_DIR=/Volumes/Alter/Daus_WGS_Paper
ALN_DIR=${HOME_DIR}/Alignments
REF_DIR=${HOME_DIR}/References
ANALYSIS_DIR=${HOME_DIR}/Analysis
SOFT_DIR=/Users/ollie
ANGSD=${SOFT_DIR}/angsd

REFERENCE=${REF_DIR}/LHISI_Scaffold_Assembly.fasta

# Make output directory
mkdir -p ${ANALYSIS_DIR}/SitewiseH
cd ${ANALYSIS_DIR}/SitewiseH

# Extract depth files
mkdir depth_files
ARCHIVE=$(find ${ANALYSIS_DIR}/SingleSampleHEstimation -name 'SingleSampleH_*.tar.gz' -type f | head -n 1)
tar -xzvf ${ARCHIVE} -C depth_files '*_nonRepeat'
# Lots of files, so it takes a bit

# Get files for coordinates of all site types
for TYPE in MODIFIER LOW MODERATE HIGH
do
  awk -v TYPE=${TYPE} 'NR >1 && $4 == TYPE {print $1"\t"$2-1"\t"$2}' ${ANALYSIS_DIR}/DeleteriousMutations/vars_classified.txt > ${TYPE}.bed
done

# Now loop over all individuals, depths, chromosomes, and site types
# Define arrays
IND_ARR=(C01210 C01211 C01213 C01215 C01217 C01219 C01220 C01222 C01223 C01224 C01225 C01227 C01230 C01231 C01232 C01234 C10223 PAUL4 VAN2)
DEPTH_ARR=(low mid)
CHROM_ARR=(CM056993.1 CM056994.1 CM056995.1 CM056996.1 CM056997.1 CM056998.1 CM056999.1 CM057000.1 CM057001.1 CM057002.1 CM057003.1 CM057004.1 CM057005.1 CM057006.1 CM057007.1 CM057008.1)
TYPE_ARR=(HIGH MODERATE LOW MODIFIER)

IND_ARR=(C01210)
DEPTH_ARR=(low)
CHROM_ARR=(CM056993.1 CM056994.1 CM056995.1 CM056996.1 CM056997.1 CM056998.1 CM056999.1 CM057000.1 CM057001.1 CM057002.1 CM057003.1 CM057004.1 CM057005.1 CM057006.1 CM057007.1 CM057008.1)
TYPE_ARR=(MODIFIER)

for IND in "${IND_ARR[@]}"; do
  for DEPTH in "${DEPTH_ARR[@]}"; do
    for TYPE in "${TYPE_ARR[@]}"; do

      #Prefix to save time
      PREFIX=${IND}_${DEPTH}_${TYPE}

      [ -f "${PREFIX}.bed" ] && rm "${PREFIX}.bed"
      [ -f "regions.bed" ] && rm "regions.bed"

      # Loop over all chroms and paste sites into bed file
      for CHROM in "${CHROM_ARR[@]}"; do

        # Get overlapping sites at chrom, depth, and type
        bedtools intersect -a ${TYPE}.bed -b depth_files/${IND}_${DEPTH}Depth_${CHROM}_nonRepeat >> ${PREFIX}.bed

        # Also use this opportunity to make a regions file
        echo "${CHROM}:1-" >> regions.bed

      done

      # If there are no sites... remove it!
      if [[ $(wc -l < "${PREFIX}.bed") -eq 0 ]]; then
        rm ${PREFIX}.bed
        continue
      fi

      # Index the sites
      ${ANGSD}/angsd sites index ${PREFIX}.bed

      # Get sample allele frequencies per site
      ${ANGSD}/angsd \
      -i ${ALN_DIR}/${IND}.bam \
      -anc ${REFERENCE} -ref ${REFERENCE} \
      -out ${PREFIX} \
      -nThreads 8  \
      -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 30 -minQ 30 \
      -doSaf 1 -GL 2 -rf regions.bed -sites ${PREFIX}.bed 2> ${PREFIX}_saf.err

      # Get single sample SFS
      ${ANGSD}/misc/realSFS \
      -P 8 -fold 1 -maxIter 5000 \
      ${PREFIX}.saf.idx > ${PREFIX}.ml 2> ${PREFIX}_ml.err

      # Collate info into table
      FILECONTENTS=$(<"${PREFIX}.ml")
      echo -e "${IND}\t${DEPTH}\t${TYPE}\t${FILECONTENTS}" >> site_hets.txt

      # Cleanup output
      rm *.saf.idx *.saf.pos.gz *.saf.gz *arg *.ml ${PREFIX}.bed*

    done ; done ; done

# Now run a slightly different analysis
# Use angsd -doHWE to estimate persite H in LHISI and LHIP
# Then, we subset this to just the specific sites of interest later
