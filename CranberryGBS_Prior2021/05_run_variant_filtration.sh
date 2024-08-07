#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="capture seq variant calling - variant filtration"
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
PROJNAME=PROJECTNAME

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/PATH/TO/PROJECT

# Directory to output variants
VARIANTDIR=$WD/variant_calling/variants/

# Path to BED file containing probe locations
PROBEBED=/project/gifvl_vaccinium/cranberryGenotyping/PATH/TO/PROBE/BEDFILE

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
VARIANTFILES=$VARIANTDIR/${PROJNAME}_variants.vcf

# bcftools stats
OUTPUTSTAT=${VARIANTFILES%".vcf.gz"}_stats.txt
bcftools stats $VARIANTFILES > $OUTPUTSTAT


OUTPUT=${VARIANTFILES%".vcf.gz"}_filtered.vcf.gz

vcftools --gzvcf $VARIANTFILES \
	--bed $PROBEBED \
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

