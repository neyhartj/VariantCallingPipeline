#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="merged rapid variants - variant filtration"
#SBATCH -p short
#SBATCH -t 02:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 8   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## Capture Seq Variant Calling Pipeline
##
## Step 2. Filter variants from freebayes
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

# Parameters for VCFtools filtration
MAXMISSING=0.2
MINDP=10



##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Get the variant files
# VARIANTFILES=$VARIANTDIR/${PROJNAME}_variants.vcf
VARIANTFILES=$VARIANTDIR/rapid_merge_all_benlear-ref_call_20240304.vcf.gz

# bcftools stats
OUTPUTSTAT=${VARIANTFILES%".vcf.gz"}_stats.txt
bcftools stats $VARIANTFILES > $OUTPUTSTAT


## Run vcftools filter
OUTPUT=${VARIANTFILES%".vcf.gz"}_filtered.vcf.gz

vcftools --gzvcf $VARIANTFILES \
	--remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
	--max-missing $MAXMISSING \
	--minDP $MINDP \
	--recode \
	--recode-INFO-all \
	--stdout | gzip -c > $OUTPUT


# bcftools stats
OUTPUTSTAT=${OUTPUT%".vcf.gz"}_stats.txt
bcftools stats $OUTPUT > $OUTPUTSTAT

