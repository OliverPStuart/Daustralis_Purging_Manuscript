#!/bin/bash
#SBATCH --job-name=rzooroh_lhisi_wgs
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=12:00:00
#SBATCH --mem=128G
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

conda activate roh_env

#######################
#### ACTUAL SCRIPT ####
#######################

cd $ANALYSIS_DIR/GenotypeLikelihoodROH

# First rescale all beagle files to pl files

for file in $(echo *beagle*)
do
	Rscript rescale_gls.R ${file}
done

# Now tar all the beagle files because they're not
# useful anymore.

tar -cvf beagle_files.tar *beagle*
rm -rf *beagle.gz

# Now combine all of the per-scaffold results into one
# file per population.

for pop in lhisi lhip all
do
	zcat ${pop}_CM*.pl.gz >> ${pop}_genome.pl
	#gzcat ${pop}_CM*.pl.gz >> ${pop}_genome.pl
done

gzip *pl

rm *_CM*

# First, we run a script to estimate the best number
# of ROH classes to model. The package, RZooROH uses
# a HMM approach, and you need to specify the number
# of ROH classes that the model will attempt to identify
# the transitions between. We use the same set of rates
# for all models, to make the models comparable by BIC.

# We do this separately for each population subdivision.
# To change the rates and classes used, modify the Rscript.

# You may want to run this on a cluster, and manually
# uncomment the line which parallelises the model estimation.

for pop in lhisi lhip all
do
  Rscript model_comparison.R ${pop}_genome.pl.gz
done

# Now we manually observe the results and plot with other scripts
