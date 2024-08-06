#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="capture seq variant calling - alignment merging"
#SBATCH -p short
#SBATCH -t 12:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=156G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## GBS variant calling pipeline
##
## Step 3. Mark and remove duplicates in the BAM files and then merge bam files
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load samtools
module load picard

#####################
## Set variables
#####################

# A name for the project
PROJNAME="cranberry_prior2021_gbs"

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/cranberryHistoricalGBS

# Directory containing alignment files
ALIGNDIR=$WD/alignment/

# Directory to output merged alignment files
MERGEDALIGNDIR=$WD/merged_alignment

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# List files for the Stevens alignment
BAMFILES=$(find $ALIGNDIR -name "*_alignment.bam")

# Iterate over the alignment files
for BAMFILE in $BAMFILES; do

	# Output file from duplicate marking
	OUTPUT=${BAMFILE%".bam"}_nodup.bam
	MARKDUPOUT=${BAMFILE%".bam"}_duplicate_metrics.txt

	# Mark and remove duplicates
	java -Xmx100G -jar picard.jar MarkDuplicates \
		--REMOVE_DUPLICATES true \
		-I $BAMFILE \
		-O $OUTPUT \
		-M $MARKDUPOUT

done


# Collect new BAM file names
NEWBAMFILES=$(find $ALIGNDIR -name "*_alignment_nodup.bam")

# Merge the bam files
# Sort on coordinates
# Subset the bam files for only those positions overlapping with the probe BED file
samtools merge -@ $SLURM_JOB_CPUS_PER_NODE -o - $NEWBAMFILES | \
	samtools sort -@ $SLURM_JOB_CPUS_PER_NODE -u - | \
	samtools view -b -o $MERGEDALIGNDIR/${PROJNAME}_alignments_merged.bam -@ $SLURM_JOB_CPUS_PER_NODE -

# Index
samtools index $MERGEDALIGNDIR/${PROJNAME}_alignments_merged.bam
