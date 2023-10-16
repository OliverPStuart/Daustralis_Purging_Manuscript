#!/bin/bash

### Repeat masking of new Dryococelus australis assembly
### Configured to run on Dayhoff

# Get shell configuration from one directory up from script
# Move to annotation directory
# RepeatModeler doesn't let you specify output directory manually

. ../shell.config
cd ${SCRATCH_DIR}/References/Annotation

# Download singularity container

singularity pull dfam-tetools-latest.sif docker://dfam/tetools:latest

# Now run repeatmodeler

singularity exec repeatmodeler.sif BuildDatabase -name LHISI_Scaffold_Assembly_Polished
singularity exec repeatmodeler.sif RepeatModeler -database LHISI_Scaffold_Assembly_Polished -threads 28

# Now use CD-hit to reduce redundancy with the previous run

conda activate cdhit_env
cat LHISI_Scaffold_Assembly_Polished-families.fa LHISI_Scaffold_Assembly.repeats_filtered.fasta > combined_repeats.fa
cd-hit-est combined_repeats.fa -c 0.99 -o combined_repeats_clustered.fa
conda deactivate
