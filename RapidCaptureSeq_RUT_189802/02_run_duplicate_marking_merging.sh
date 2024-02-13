#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="capture seq variant calling - alignment merging"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=156G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## Capture Seq Variant Calling Pipeline
##
## Step 2. Mark and remove duplicates in the BAM files and then merge bam files
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
PROJNAME=PROJECTNAME

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/PATH/TO/PROJECT

# Directory containing alignment files
ALIGNDIR=$WD/variant_calling/alignment/

# Directory to output merged alignment files
MERGEDALIGNDIR=$WD/variant_calling/merged_alignment

# Path to BED file containing probe locations
PROBEBED=/project/gifvl_vaccinium/cranberryGenotyping/PATH/TO/PROBE/BEDFILE

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
	java -Xmx100G -jar /software/el9/apps/picard/3.0.0/picard.jar MarkDuplicates \
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
	samtools view -b -o $MERGEDALIGNDIR/${PROJNAME}_alignments_merged.bam -@ $SLURM_JOB_CPUS_PER_NODE -L $PROBEBED -

# Index
samtools index $MERGEDALIGNDIR/${PROJNAME}_alignments_merged.bam
