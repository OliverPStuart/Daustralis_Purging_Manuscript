#!/bin/bash

# Script will not source config file below
# Already checked paths and whether script knows file is there
# Not an alias issue
# Not any other recognisable issue

. ../../shell.config

# Will probably rewrite this as a snakefile!

# Whole operation took ~ 6 minutes on smallest vcf (last chromosome) on laptop WSL

# Define vcf

cd $ANALYSIS_DIR/PAUL4
VCF=PAUL4_Raw_CM057008.1.vcf
VCF_BASENAME=$(echo $VCF | sed 's/\.vcf//g' | sed s'/_Raw//g')

# Filter by depth first
# This will significantly reduce the burden of later sections
# Permissive filter

vcftools --vcf $VCF \
--min-meanDP 15 --max-meanDP 25 \
--remove-indels \
--recode --recode-INFO-all --out depthfilt_$VCF_BASENAME

# Split into header, variant, and invariant

grep "#" depthfilt_${VCF_BASENAME}.recode.vcf > header
grep -v "#" depthfilt_${VCF_BASENAME}.recode.vcf > tailer
rm depthfilt_${VCF_BASENAME}.recode.vcf

grep "0/0\|1/1" tailer > tmp1

grep "0/1\|1/0" tailer > tmp2
cat header tmp2 > variant_depthfilt_$VCF_BASENAME.vcf
rm tmp2

rm tailer *log

# Filter variant VCF, easy part

vcffilter -f 'AB > 0.2 & AB < 0.8 & QUAL / DP > 0.25' variant_depthfilt_$VCF_BASENAME.vcf | grep -v "#" >  tmp2

# These may have been true homozygotes
# The imbalanced filtering of het/hom sites due to information imbalance worries me
# So I'm recording how much information is lost when hets are filtered
# Perhaps I can scale down the number of homozygous sites by a proportional amount to account for this
# Or perhaps I can calculate metrics assuming these removed sites were true homozygotes

NUM_REMAINING=$(cat tmp2 | wc -l)
NUM_TOTAL=$(grep -v "#" variant_depthfilt_$VCF_BASENAME.vcf | wc -l)
echo "${VCF_BASENAME}\t${NUM_REMAINING}\t${NUM_TOTAL}" >> het_filtering.log

# Now cat together, sort

cat tmp1 tmp2 | sort -k2,2n > tmp3
rm tmp1 tmp2

cat header tmp3 > FILTERED_$VCF_BASENAME.vcf
rm tmp3 header variant_depthfilt_$VCF_BASENAME.vcf

