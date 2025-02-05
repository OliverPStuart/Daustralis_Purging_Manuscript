# Purging and genetic diversity in *Dryococelus australis*
All scripts used to analyse genetic diversity in the Lord Howe Island stick insect *Dryococelus australis*. All reads used for this are available under BioProject PRJNA1124435, totally 24 individuals and 25 sequence datasets. Below is a rough roadmap of the separate analyses run. Software requirements are not described.

## 1 High-coverage data
Analysis of one wild male individual (Paul) sequenced to high (~20x) coverage.
- `Scripts/SLURM/reference_index_slurm.sh` indexes the reference genome, required for future steps.
- `Scripts/Snakefiles/PAUL4.smk` performs all bioinformatic processing of Paul's sequence data, alignment, indexing, variant calling, depth analysis, heterozygosity, and ROH with bcftools. Parallelised across chromosomes or, for variant calling, over many chunks of a defined number per chromosome.
- `Scripts/Shell/collate_depth.sh` combines per-sample depth estimates into a single table.
- `Scripts/Shell/paul4_gc_content.sh` estimates GC content for analysis.
- `Scripts/RScripts/paul4_analysis.R` analyses data generated in the previous step, generating figures.

## 2 Low-coverage data
Analysis of 24 individuals from wild, captive, and hybrid individuals using low-coverage (~2-6x) sequencing. Paul is included among the wild individuals but with different sequence data.

### 2.1 Alignment and QC
- `Scripts/Snakefiles/AlignmentWGS.smk` combines read files per sample and aligns/indexes, also estimating depth.
- `Scripts/RScripts/wgs_depth_analysis.R` analyses depth data and generates figures.

### 2.2 Variant identification
- `Scripts/Snakefiles/AlleleFrequencies.smk` estimates allele frequencies across all low-coverage individuals, then takes every identified site and estimates group specific frequencies.
- `Scripts/Shell/ngsparalog_analysis.sh` takes allele frequency estimates and runs ngsParalog to identify regions of putative mismapping, used in later steps for filtering.

### 2.3 Population structure
- `Scripts/SLURM/lcwgs_pcangsd_slurm.sh` estimates sample covariance matrices for all samples accounting for masked regions.
- `Scripts/RScripts/pcangsd_analysis.R` analyses covariance matrices and generates figures.

### 2.4 Individual heterzygosity
- `Scripts/Snakefiles/SingleSampleHEstimation.smk` estimates per sample heterozygosity accounting for masked regions and inter-individual depth variation.
- `Scripts/RScripts/singlesampleH_analysis.R` analyses heterozygosity values and generates figures.

### 2.5 Runs-of-homozygosity
- `Scripts/Snakefiles/WGSGenotypeLikelihoods.smk` estimates genotype likelihoods at variable positions for all individuals accounting for masked regions.
- `Scripts/lcwgs_rzooroh_slurm.sh` estimates the locations of ROH for all individuals. Refers to the following scripts:
  - `Scripts/RScripts/rescale_gls_for_rzooroh.R` rescales genotype likelihood files from the first step to Phred scale.
  - `Scripts/RScripts/rzooroh_model_comparison.R` runs the model estimation procedure, producing summary figures of BIC/AIC/LLs for different parameter combinations.
- `Scripts/RScripts/genotype_likelihood_roh_analysis.R` analyses ROH results and generates figures.
- `Scripts/RScripts/population_base_bias_in_f_estimate.R` analyses ROH results from different base population models and generates figures.

### 2.6 Deleterious mutations
- `Scripts/Shell/build_snpeff_database.sh` builds a SnpEff database from the *D. australis references genome*
- `Scripts/classify_var.sh` takes positions and non-reference alleles from allele frequency estimation and estimates their effects on protein function using the SnpEff database. Also classifies positions by how many individuals have that position contained within ROH.
- `Scripts/RScripts/deleterious_mutation_analysis.R` analyses the deleterious mutation data and generates figures.
