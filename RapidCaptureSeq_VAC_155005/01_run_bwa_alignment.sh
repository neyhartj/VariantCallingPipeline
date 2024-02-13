#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="VAC_155005 - capture seq variant calling - alignment"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## Capture Seq Variant Calling Pipeline
##
## Step 1. Alignment to reference genomes using the BWA aligner
##

# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load bwa
module load samtools

#####################
## Set variables
#####################

# A name for the project
PROJNAME=VAC_155005

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/RAPiD_Cranberry_15K/Data/2022/

# Path to directory containing the FASTQ files
FASTQDIR=$WD/fastq_files/

# Path to TXT file containing sample names to run (one per line)
SAMPLEFILE=$WD/VAC_155005_sample_file.txt

# Prefix of the indexed reference genome
REFPREFIX=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_macrocarpon_Stevens_v1.fasta

# Directory to output alignment files
ALIGNDIR=$WD/variant_calling/alignment/

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Use the sample file to create a vector of sample names
SAMPLENAMES=$(cut -d \t -f 1 $SAMPLEFILE)

# Iterate over the sample names
for SAMPLE in $SAMPLENAMES; do
  # Create a RG tag
  RG="@RG\tID:$SAMPLE\tSM:$SAMPLE"

  # Find the FASTQ files in the input directory that match the sample name
  SAMPLEFASTQS=$(find $FASTQDIR -name "*$SAMPLE*")

  ## ALIGNMENT TO STEVENS
  # Create the output SAM file name
  OUTPUT=$ALIGNDIR/${SAMPLE}_alignment.sam
  OUTPUTBAM=$ALIGNDIR/${SAMPLE}_alignment.bam
  
  # Run the alignment in a pipeline
  bwa mem -t $NTHREADS -R $RG $REFPREFIX $SAMPLEFASTQS | \
  samtools fixmate -u -m - - | \
  samtools sort -b -@ $NTHREADS -o $OUTPUTBAM -

done
