#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=08:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="stacks_process_radtags"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## STACKS Process Radtags pipeline
## 
## This script will run the process_radtags function from stacks to demultiplex,
## cleaning, and quality control
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load stacks

## Set variables

# Working directory
WD=/project/gifvl_vaccinium/cranberryGenotyping/cranberryHistoricalGBS

# Name of input directory
INPUT=$WD/input
# Name of input directory with FASTQ
FASTQDIR=$INPUT/fastq_files
# Name of the input directory with the barcode files
BCODEDIR=$INPUT/stacks_barcode_files
# Name of the output directory
OUTPUT=$WD/demultiplexed_fastq_files

# Number of threads available
NTHREADS=$SLURM_JOB_CPUS_PER_NODE


#########################################################
# Run the pipeline

# Change working directory
cd $WD

# List the barcode files
bcodefiles=$(find $BCODEDIR -name *stacks_barcodes.txt)

# Iterate over the barcode files
for file in $bcodefiles; do
    
    # Identify the flowcell lane for this file
    flowcell=$(basename $file | sed 's/_stacks_barcodes.txt//g')
    
    # Make a folder in the output dir for this flowcell
    outputdir=$OUTPUT/$flowcell
    mkdir $outputdir
    
    # Find that fastq file in the input directory
    fastqfile=$(find $FASTQDIR -name ${flowcell}_fastq.gz)
    
    # Run process radtags
    process_radtags -f $fastqfile -b $file -e 'ecoT22I' -o $outputdir --threads $NTHREADS --clean --quality --rescue
    
done
    
    