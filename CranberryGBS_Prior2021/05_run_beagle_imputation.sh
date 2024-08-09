#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=04:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=40G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="beagle_imputation"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## Beagle imputation testing
## 
## This script will run Beagle imputation
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load beagle-geno
module load vcftools

#####################
## Set variables
#####################

# A name for the project
PROJNAME="cranberry_prior2021_gbs"

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/cranberryHistoricalGBS

# Directory to output variants
VARIANTDIR=$WD/variants/
# Name of input vcf
VCFIN=$VARIANTDIR/cranberry_prior2021_gbs_variants_filtered_imputation.vcf.gz

# Output directory
OUTPUT=$WD/imputation/beagle_imputation


# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE



##############################
## DO NOT EDIT BELOW
##############################


# Change working directory
cd $WD 

## Run beagle

# Run the imputation software
java -Xmx32g -jar /software/7/apps/beagle-geno/5.3/beagle.08Feb22.fa4.jar \
gt=$VCFIN \
out=$OUTPUT/${PROJNAME}_imputed_snps \
ne=25 burnin=10 iterations=30 nthreads=$NTHREADS

