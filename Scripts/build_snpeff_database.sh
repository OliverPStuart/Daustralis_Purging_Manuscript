### Building snpEff database for Dryococelus australis

# Download snpeff
cd ~
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

# Download gffread
cd ~
git clone https://github.com/gpertea/gffread
cd gffread
make release
cd ..

# Add D. australis to snpEff config file
cd snpEff
echo "# Dryococelus australis Daus2.0" >> snpEff.config
echo "Daus2.0.genome: Dryococelus australis" >> snpEff.config
echo "" >> snpEff.config

# Make data folder and download genome
mkdir data/Daus2.0
cd data/Daus2.0

conda activate ncbi_download_env
datasets download genome accession GCA_029891345.1 --include genome
conda deactivate

unzip ncbi_dataset.zip
mv ncbi_dataset/data/GCA_029891345.1/GCA_029891345.1_Daus_2.0_genomic.fna sequences.fa
rm -rf ncbi_dataset*


# Use original annotation, but have to rename the scaffolds in there
pip install osfclient
osf -p 9v487 fetch osfstorage/LHISI_Scaffold_Assembly.annotation.gff.gz
gunzip LHISI_Scaffold_Assembly.annotation.gff.gz

echo -e "CM056993.1\tScaffold_1\nCM056994.1\tScaffold_2\nCM056995.1\tScaffold_3\nCM056996.1\tScaffold_5\nCM056997.1\tScaffold_6\nCM056998.1\tScaffold_7\nCM056999.1\tScaffold_8\nCM057000.1\tScaffold_9\nCM057001.1\tScaffold_10\nCM057002.1\tScaffold_11\nCM057003.1\tScaffold_12\nCM057004.1\tScaffold_13\nCM057005.1\tScaffold_14\nCM057006.1\tScaffold_15\nCM057007.1\tScaffold_16\nCM057008.1\tScaffold_17\nCM057009.1\tScaffold_4\n" > scaffold_rename_table

awk 'NR==FNR{a[$2]=$1} NR>FNR{$1=a[$1];print}' OFS="\t" \
scaffold_rename_table \
LHISI_Scaffold_Assembly.annotation.gff > genes.gff

rm LHISI_Scaffold_Assembly.annotation.gff

cut -f9- genes.gff | sed -e "s/\t/ /g" > col9_reformat
cut -f1-8 genes.gff | paste - col9_reformat > tmp
mv tmp genes.gff

# Get cds and proteins from the  annotation
~/gffread/gffread -g sequences.fa -x cds.fa genes.gff
~/gffread/gffread -g sequences.fa -y protein.fa genes.gff

# Build database
cd ~/snpEff
java -jar snpEff.jar build -gff3 -v Daus2.0

# Should take 1-2 minutes
