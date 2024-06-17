### Script to analyses GC content of PAUL4's genome

# Set directories

HOME_DIR=/mnt/data/dayhoff/home/scratch/groups/mikheyev/LHISI
#HOME_DIR=/Volumes/Alter/Daus_WGS_Paper
REF_DIR=${HOME_DIR}/References
WORKING_DIR=${HOME_DIR}/Analysis/PAUL4
BCFTOOLS=/mnt/data/dayhoff/home/u6905905/bcftools/bcftools

cd ${WORKING_DIR}

conda activate paul4_env

# Use PAUL4 1/1 SNP calls to modify genome fasta to get PAUL4 consensus sequence

grep "#" PAUL4_Hom_CM056993.1.vcf > head

for i in *Hom*
do
grep "1/1" ${i} >> tail
done

cat head tail > PAUL4_Genome_HomAlt.vcf
rm head tail

bgzip PAUL4_Genome_HomAlt.vcf
tabix PAUL4_Genome_HomAlt.vcf.gz
cat ${REF_DIR}/LHISI_Scaffold_Assembly.fasta | \
${BCFTOOLS} consensus PAUL4_Genome_HomAlt.vcf.gz > PAUL4_Consensus.fasta

# Generate random 100 bp intervals on callable genome
# This is genome minus repeats
# Remove any intervals

awk '{print $1 "\t" $3+1}' ${REF_DIR}/autosome_regions.bed > lhisi.genome

bedtools random -g lhisi.genome > random_regions.bed

bedtools subtract \
-a random_regions.bed \
-b ${REF_DIR}/Annotation/repeats.bed > masked_random_regions.bed

# Analyse GC content of these random regions

echo -e "chr\tlength\tA\tC\tG\tT" > paul4_gc.txt
seqtk comp -r masked_random_regions.bed PAUL4_Consensus.fasta | \
awk '{print $1"\t"$3-$2"\t"$4"\t"$5"\t"$6"\t"$7}' >> paul4_gc.txt

# Analyse in R
