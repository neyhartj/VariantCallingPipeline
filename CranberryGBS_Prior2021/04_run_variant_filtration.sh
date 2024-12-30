#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="cranberry GBS variant calling - variant filtration"
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
PROJNAME="cranberry_prior2021_gbs"
# REFNAME="stevensv1"
REFNAME="benlearv1"

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/cranberryHistoricalGBS

# Directory to output variants
VARIANTDIR=$WD/variants/

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE

# Filtering parameters for production SNPs
MAXMISSING=0.50
MINDP=5
MINMAC=1 # Remove monomorphic SNPs
MINQ=10 # Minimum QUAL score


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Get the variant files
VARIANTFILES=$VARIANTDIR/${PROJNAME}_${REFNAME}_variants.vcf.gz

# Run stats
OUTPUTSTAT=${VARIANTFILES%".vcf.gz"}_stats.txt
bcftools stats $VARIANTFILES > $OUTPUTSTAT

### Filter for production SNPS ##
OUTPUT=${VARIANTFILES%".vcf.gz"}_filtered

vcftools --gzvcf $VARIANTFILES \
	--remove-indels \
	--min-alleles 2 \
	--max-alleles 2 \
	--max-missing $MAXMISSING \
	--mac $MINMAC \
	--recode \
	--recode-INFO-all \
	--out $OUTPUT
	
# Additional filter with bcftools for allele depth and quality
VARIANTFILE1=${OUTPUT}.recode.vcf

VARIANTFILE2=${OUTPUT}.vcf
bcftools view -i "AD[GT] >= $MINDP && F_MISSING <= $MAXMISSING" $VARIANTFILE1 > $VARIANTFILE2

# Zip
gzip -f $VARIANTFILE2

# Run stats on the output
OUTPUTSTAT=${VARIANTFILE2%".vcf.gz"}_stats.txt
bcftools stats $VARIANTFILE2 > $OUTPUTSTAT



# ### Filter for haplotype library construction ##
# 
# OUTPUT=${VARIANTFILES%".vcf.gz"}_filtered_haplotype_library
# 
# vcftools --gzvcf $VARIANTFILES \
# 	--remove-indels \
# 	--min-alleles 2 \
# 	--max-alleles 2 \
# 	--max-missing 0.80 \
# 	--mac 5 \
# 	--minDP 5 \
# 	--keep $VARIANTDIR/cranberry_gbs_imputation_individuals.txt \
# 	--chr chr1 \
# 	--chr chr2 \
# 	--chr chr3 \
# 	--chr chr4 \
# 	--chr chr5 \
# 	--chr chr6 \
# 	--chr chr7 \
# 	--chr chr8 \
# 	--chr chr9 \
# 	--chr chr10 \
# 	--chr chr11 \
# 	--chr chr12 \
# 	--recode \
# 	--recode-INFO-all \
# 	--out $OUTPUT
# 	
# 
# # Create a BED file based on this VCF file
# cat $OUTPUT.recode.vcf | grep -v "^#" | awk 'NR > 1 { print $1"\t"$2"\t"$2 }' - > ${OUTPUT}_snps.bed
# 
# BEDFILE=${OUTPUT}_snps.bed
# 
# ### Filter for imputation ##
# 
# OUTPUT=${VARIANTFILES%".vcf.gz"}_filtered_imputation
# 
# vcftools --gzvcf $VARIANTFILES \
# 	--remove-indels \
# 	--min-alleles 2 \
# 	--max-alleles 2 \
# 	--mac 5 \
# 	--minDP 5 \
# 	--bed $BEDFILE \
# 	--recode \
# 	--recode-INFO-all \
# 	--out $OUTPUT
# 
