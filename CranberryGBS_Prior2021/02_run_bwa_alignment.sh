#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --job-name="GBS variant calling - alignment"
#SBATCH -p short
#SBATCH -t 24:00:00   # walltime limit (HH:MM:SS)
#SBATCH -N 1   # number of nodes
#SBATCH -n 32   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=8G   # maximum memory per node
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN,END,FAIL


##
## GBS variant calling pipeline
##
## Step 2. Alignment to reference genomes using the BWA aligner
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
PROJNAME="cranberry_prior2021_gbs"

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/cranberryHistoricalGBS

# Path to directory containing the FASTQ files
FASTQDIR=$WD/demultiplexed_fastq_files

# Prefix of the indexed reference genome
REFPREFIX=/project/gifvl_vaccinium/cranberryGenotyping/genome_assemblies/Vaccinium_macrocarpon_Stevens_v1.fasta

# Directory to output alignment files
ALIGNDIR=$WD/alignment/

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


##############################
## DO NOT EDIT BELOW
##############################


## Run the pipeline

# Change working directory
cd $WD

# Make a directory for the merged fastq files
FASTQMERGEDIR=$WD/merged_fastq_files
if [ ! -d $FASTQMERGEDIR ]; then
  mkdir $FASTQMERGEDIR
fi

# # Find fastq files
# FASTQFILES=$(find $FASTQDIR -name "*.fq.gz")
# # Get unique sample names
# FASTQUNIQUE=$(for file in $FASTQFILES; do basename $file; done | sort -u)
# 
# # Iterate over the unique sample names and merge duplicates
# for sample in $FASTQUNIQUE; do
#   # Search for the sample in FASTQFILES
#   sample_fastqs=($(find $FASTQDIR -name $sample))
#   outname=$FASTQMERGEDIR/$sample
#   # If more than one, merge
#   if [ "${#sample_fastqs[@]}" -gt 1 ]; then
#     cat ${sample_fastqs[@]} > $outname
#   else
#     cp ${sample_fastqs[@]} $outname
#   fi
# done
  

## Get a new list of the fastq files
FASTQFILES=$(find $FASTQMERGEDIR -name "*.fq.gz")

# Iterate over those files and align to the Stevens reference genome
for fastqfile in $FASTQFILES; do
  SAMPLE=$(basename $fastqfile | sed 's/.fq.gz//g')
  # Create a RG tag
  RG="@RG\tID:$SAMPLE\tSM:$SAMPLE"
  ## ALIGNMENT TO STEVENS
  # Create the output SAM file name
  OUTPUT=$ALIGNDIR/${SAMPLE}_stevensv1_alignment.sam
  OUTPUTBAM=$ALIGNDIR/${SAMPLE}_stevensv1_alignment.bam
  
  # Run the alignment in a pipeline
  bwa mem -t $NTHREADS -R $RG $REFPREFIX $fastqfile | \
  samtools fixmate -u -m - - | \
  samtools sort -O bam -@ $NTHREADS -o $OUTPUTBAM -
    
done

