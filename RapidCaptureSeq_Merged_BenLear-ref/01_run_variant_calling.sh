#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="merged rapid variants - variant calling"
#SBATCH -p short
#SBATCH -t 06:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 8   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## Capture Seq Variant Calling Pipeline
##
## Step 4. Filter variants from freebayes
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load vcftools
module load bcftools

#####################
## Set variables
#####################

# A name for the project
PROJNAME=Rapid_Merged_BenLear-ref

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/RAPiD_Cranberry_15K/Data/2023/

# Directory to output variants
VARIANTDIR=$WD/variants/

# Path to BED file containing probe locations
PROBEBED=""

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Get the variant files
# VARIANTFILES=$VARIANTDIR/${PROJNAME}_variants.vcf
VARIANTFILES=$VARIANTDIR/rapid_merge_all_benlear-ref_20240223.vcf.gz

# Run bcftools call
OUTPUT=${VARIANTFILES%".vcf.gz"}_call.vcf.gz
bcftools call $VARIANTFILES -o $OUTPUT -O z -mv --threads $NTHREADS


