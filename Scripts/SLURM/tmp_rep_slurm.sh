#!/bin/bash
#SBATCH --job-name=pcangsd_lhisi_wgs_rep
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=36:00:00
#SBATCH --mem=36G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

###################
#### ENV SETUP ####
###################

### Source .bashrc for conda environments
#source ${HOME_DIR}/.bashrc
eval "$(conda shell.bash hook)"

### Evironment
cd /mnt/data/dayhoff/home/u6905905/Daus_WGS_Paper
. shell.config

SOFTWARE_DIR=/mnt/data/dayhoff/home/u6905905
PCANGSD=${SOFTWARE_DIR}/pcangsd/pcangsd/pcangsd.py
ANGSD=${SOFTWARE_DIR}/angsd/angsd

REFERENCE=${REF_DIR}/LHISI_Scaffold_Assembly.fasta

conda activate pcangsd_env

#######################
#### ACTUAL SCRIPT ####
#######################

cd ${SCRATCH_DIR}/Analysis
#mkdir -p PCangsd
cd PCangsd

# Make a file containing the names of all individuals to be included
# Remove samples that didn't run well

#find $ALN_DIR | grep bam$ | grep -v "C01216\|C01233\|C01218\|C01226\|C10133\|PAUL4_WGS" > all_bam_list.txt

# Now loop over the estimated allele frequency files in the AlleleFrequency
# analysis directory. We want to get lists of sites in the genome estimated as
# variant previously

# First get file of all sites

#for scaffold in $(cat $REF_DIR/autosome_names)
#	do
#	zcat ${SCRATCH_DIR}/Analysis/AlleleFrequencies/${scaffold}_all.mafs.gz | \
#	awk 'NR > 1 {print $1"\t"$2-1"\t"$2}' >> tmp1.bed
#done

# In every scaffold, we print the minor allele frequency table, then we take
# away anything that has a maf within 0.02 of 1/2N. ANGSD will sometimes provide
# interesting maf estimates due to the probabilistic nature of the software e.g.
# for a sample of four chromosomes, an estimated maf might be 0.2578 instead of
# 0.25. For all samples, 1/2N is 0.026, hence the 0.02 threshold.

#for scaffold in $(cat $REF_DIR/autosome_names)
#	do
#	zcat ${SCRATCH_DIR}/Analysis/AlleleFrequencies/${scaffold}_all.mafs.gz | \
#	awk 'NR > 1' | awk '($6 - 1/(2*$8)) > 0.02 || ($6 - 1/(2*$8)) < -0.02 {print $1"\t"$2-1"\t"$2}' >> tmp2.bed
#done

# When estimating allele frequencies, we only masked repetitive regions. We used
# the allele frequency site outputs to find regions of problematic mapping with
# ngsParalog. There is now a file in ${REF_DIR} called problematic_regions.bed.
# We need to use this to masked the variants.

#bedtools subtract -a tmp1.bed -b ${REF_DIR}/problematic_regions.bed | \
#awk '{print $1"\t"$3}' > sites
#bedtools subtract -a tmp2.bed -b ${REF_DIR}/problematic_regions.bed | \
#awk '{print $1"\t"$3}' > sites_maf

#rm tmp*

# Now format as a vcf
#echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > dummy.vcf
#awk '{print $1"\t"$2"\t"$1"_"$2"\tA\tT\t100.00\tPASS\tNA"}' sites >> dummy.vcf

#echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > dummy_maf.vcf
#awk '{print $1"\t"$2"\t"$1"_"$2"\tA\tT\t100.00\tPASS\tNA"}' sites_maf >> dummy_maf.vcf

# Now we can thin this to get pseudo-independence
#vcftools --vcf dummy.vcf --thin 50000 --recode --recode-INFO-all --out thin
#vcftools --vcf dummy_maf.vcf --thin 50000 --recode --recode-INFO-all --out thin_maf

# Print as sites again
#awk 'NR > 1 {print $1,$2}' thin.recode.vcf > sites_thin
#awk 'NR > 1 {print $1,$2}' thin_maf.recode.vcf > sites_maf_thin
#rm *vcf *log

# Now index these with angsd

#${ANGSD} sites index sites_thin
#${ANGSD} sites index sites_maf_thin

#sleep 10s

#touch sites_thin.*
#touch sites_maf_thin.*

# Now use angsd to get a genotype likelihood file for these sites

# Test line for small file
#${ANGSD} -b all_bam_list.txt -ref ${REFERENCE} -out all_thin \
#-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 30 -minQ 30 \
#-minInd 9 -setMinDepthInd 1 -setMaxDepthInd 8 \
#-doCounts 1 -GL 2 -doGlf 2 -doMajorMinor 4 -nThreads 8 \
#-r CM057008.1:1-1000000

#${ANGSD} -b all_bam_list.txt -ref ${REFERENCE} -out all_thin \
#-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 30 -minQ 30 \
#-minInd 9 -setMinDepthInd 1 -setMaxDepthInd 8 \
#-doCounts 1 -GL 2 -doGlf 2 -doMajorMinor 4 -nThreads 8 \
#-sites sites_thin

${ANGSD} -b all_bam_list.txt -ref ${REFERENCE} -out all_thin_maf \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 30 -minQ 30 \
-minInd 9 -setMinDepthInd 1 -setMaxDepthInd 8 \
-doCounts 1 -GL 2 -doGlf 2 -doMajorMinor 4 -nThreads 8 \
-sites sites_maf_thin

# Now run pcangsd with different maf cutoffs and different thinning values

### Error is here... command not found
### Probably to do with python library links etc

#${PCANGSD} -b all_thin.beagle.gz -o pcangsd_all_thin --inbreed_samples
#${PCANGSD} -b all_thin.beagle.gz -o pcangsd_all_thin_05 --inbreed_samples --maf 0.05
#${PCANGSD} -b all_thin.beagle.gz -o pcangsd_all_thin_10 --inbreed_samples --maf 0.10
#${PCANGSD} -b all_thin.beagle.gz -o pcangsd_all_thin_15 --inbreed_samples --maf 0.15

#${PCANGSD} -b all_thin_maf.beagle.gz -o pcangsd_all_thin_maf --inbreed_samples
#${PCANGSD} -b all_thin_maf.beagle.gz -o pcangsd_all_thin_maf_05 --inbreed_samples --maf 0.05
#${PCANGSD} -b all_thin_maf.beagle.gz -o pcangsd_all_thin_maf_10 --inbreed_samples --maf 0.10
#${PCANGSD} -b all_thin_maf.beagle.gz -o pcangsd_all_thin_maf_15 --inbreed_samples --maf 0.15
