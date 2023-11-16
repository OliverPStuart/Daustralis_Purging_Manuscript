#!/bin/bash
#SBATCH --job-name=repeat_modeler_LHISI
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --mail-user=u6905905@anu.edu.au
#SBATCH --mail-type=FAIL,END
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

### Source .bashrc for conda environments
#source ${HOME_DIR}/.bashrc
eval "$(conda shell.bash hook)"

### Evironment
cd /mnt/data/dayhoff/home/u6905905/Daus_WGS_Paper
. shell.config

cd ${REF_DIR}

# Download singularity container
#singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest

# Now run repeatmodeler

#singularity exec -H $(pwd) dfam-tetools-latest.sif BuildDatabase -name LHISI_Scaffold_Assembly_Polished LHISI_Scaffold_Assembly_Polished.fasta
#singularity exec -H $(pwd) dfam-tetools-latest.sif RepeatModeler -database LHISI_Scaffold_Assembly_Polished -threads 28

# Split de novo ID repeats into known and unknown
#grep ">" RM_3767039.TueOct241448172023/consensi.fa.classified | grep "#Unknown" | sed 's/^>//g' > consensi.fa.classified.unknown
#grep ">" RM_3767039.TueOct241448172023/consensi.fa.classified | grep -v "#Unknown" | sed 's/^>//g' > consensi.fa.classified.known

# Now mask reference with four rounds of RepeatMasker

REFERENCE=${REF_DIR}/LHISI_Scaffold_Assembly_Polished.fasta

#mkdir -p logs 01_simple_out 02_arthropoda_out 03_known_out 04_unknown_out

conda activate repeatmasker_env

# round 1: annotate/mask simple repeats
#RepeatMasker -pa 16 -a -e ncbi -dir 01_simple_out -noint -xsmall LHISI_Scaffold_Assembly_Polished.fasta 2>&1 | tee logs/01_simplemask.log
# round 1: rename outputs
#rename 's/fasta/simple_mask/g' 01_simple_out/LHISI_Scaffold_Assembly_Polished*
#rename 's/\.masked/\.masked\.fasta/g' 01_simple_out/LHISI_Scaffold_Assembly_Polished*

# round 2: annotate/mask Arthropoda elements sourced from Repbase using output from 1st round of RepeatMasker
#RepeatMasker -pa 16 -a -e ncbi -dir 02_arthropoda_out -nolow \
#-species arthropoda 01_simple_out/LHISI_Scaffold_Assembly_Polished.simple_mask.masked.fasta 2>&1 | tee logs/02_arthropodamask.log
# round 2: rename outputs
#rename 's/simple_mask\.masked\.fasta/arthropoda_mask/g' 02_arthropoda_out/LHISI_Scaffold_Assembly_Polished*
#rename 's/\.masked/\.masked\.fasta/g' 02_arthropoda_out/LHISI_Scaffold_Assembly_Polished*

# round 3: annotate/mask known elements sourced from species-specific de novo repeat l'ibrary
RepeatMasker -pa 16 -a -e ncbi -dir 03_known_out -nolow \
-lib consensi.fa.classified.known \
02_arthropoda_out/LHISI_Scaffold_Assembly_Polished.arthropoda_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log
# round 3: rename outputs
rename 's/arthropoda_mask\.masked\.fasta/known_mask/g' 03_known_out/LHISI_Scaffold_Assembly_Polished*
rename 's/\.masked/\.masked\.fasta/g' 03_known_out/LHISI_Scaffold_Assembly_Polished*

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library
RepeatMasker -pa 16 -a -e ncbi -dir 04_unknown_out -nolow \
-lib consensi.fa.classified.unknown \
03_known_out/LHISI_Scaffold_Assembly_Polished.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log
# round 4: rename outputs
rename 's/known_mask\.masked\.fasta/unknown_mask/g' 04_unknown_out/LHISI_Scaffold_Assembly_Polished*
rename 's/\.masked/\.masked\.fasta/g' 04_unknown_out/LHISI_Scaffold_Assembly_Polished*

# create directory for full results
mkdir -p 05_full_out

# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/LHISI_Scaffold_Assembly_Polished.simple_mask.cat.gz \
02_arthropoda_out/LHISI_Scaffold_Assembly_Polished.arthropoda_mask.cat.gz \
03_known_out/LHISI_Scaffold_Assembly_Polished.known_mask.cat.gz \
04_unknown_out/LHISI_Scaffold_Assembly_Polished.unknown_mask.cat.gz \
> 05_full_out/LHISI_Scaffold_Assembly_Polished.full_mask.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/LHISI_Scaffold_Assembly_Polished.simple_mask.out \
<(cat 02_arthropoda_out/LHISI_Scaffold_Assembly_Polished.arthropoda_mask.out | tail -n +4) \
<(cat 03_known_out/LHISI_Scaffold_Assembly_Polished.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/LHISI_Scaffold_Assembly_Polished.unknown_mask.out | tail -n +4) \
> 05_full_out/LHISI_Scaffold_Assembly_Polished.full_mask.out

# copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/LHISI_Scaffold_Assembly_Polished.simple_mask.out > 05_full_out/LHISI_Scaffold_Assembly_Polished.simple_mask.out

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat 02_arthropoda_out/LHISI_Scaffold_Assembly_Polished.arthropoda_mask.out \
<(cat 03_known_out/LHISI_Scaffold_Assembly_Polished.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/LHISI_Scaffold_Assembly_Polished.unknown_mask.out | tail -n +4) \
> 05_full_out/LHISI_Scaffold_Assembly_Polished.complex_mask.out

# combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_outLHISI_Scaffold_Assembly_Polished.simple_mask.align \
02_arthropoda_out/LHISI_Scaffold_Assembly_Polished.arthropoda_mask.align \
03_known_out/LHISI_Scaffold_Assembly_Polished.known_mask.align \
04_unknown_out/LHISI_Scaffold_Assembly_Polished.unknown_mask.align \
> 05_full_out/LHISI_Scaffold_Assembly_Polished.full_mask.align

# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species arthropoda 05_full_out/LHISI_Scaffold_Assembly_Polished.full_mask.cat.gz 2>&1 | tee logs/05_fullmask.log

#rmOutToGFF3custom -o 05_full_out/LHISI_Scaffold_Assembly_Polished.full_mask.out > 05_full_out/LHISI_Scaffold_Assembly_Polished.full_mask.gff3
#rmOutToGFF3custom -o 05_full_out/LHISI_Scaffold_Assembly_Polished.simple_mask.out > 05_full_out/LHISI_Scaffold_Assembly_Polished.simple_mask.gff3
#rmOutToGFF3custom -o 05_full_out/LHISI_Scaffold_Assembly_Polished.complex_mask.out > 05_full_out/LHISI_Scaffold_Assembly_Polished.complex_mask.gff3
