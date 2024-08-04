#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=128G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="stacks_process_radtags_blueberry_GBS"
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
module load parallel

## Set variables

# Working directory
WD=/project/gifvl_vaccinium/blueberryGenotyping/blueberryEAA

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
bcodefiles=($(find $BCODEDIR -name *stacks_barcodes.txt))

# Get the flowcell_lane names from the barcode files and make directories
for file in $bcodefiles; do
  # Identify the flowcell lane for this file
  flowcell=$(basename $file | sed 's/_stacks_barcodes.txt//g')
  
  # Make a folder in the output dir for this flowcell
  # Only do this if the directory does not already exit
  outputdir=$OUTPUT/$flowcell
  if [ ! -d $outputdir ]; then
    mkdir $outputdir
  fi
  
done
  

# Write a function that takes the bcodefiles and runs the demultiplenxing function
run_demultiplex() {
  file=$1
  fastqdir=$2
  output=$3
  # Identify the flowcell lane for this file
  flowcell=$(basename $file | sed 's/_stacks_barcodes.txt//g')
  
  # Make a folder in the output dir for this flowcell
  outputdir=$output/$flowcell

  # Find that fastq file in the input directory
  fastqfile=$(find $fastqdir -name ${flowcell}_fastq.gz)
  
  # Run process radtags
  process_radtags -f $fastqfile -b $file -e 'apeKI' -o $outputdir --threads 1 --clean --quality --rescue
}

# Export the function
export -f run_demultiplex

# Run the function in parallel
parallel -j $NTHREADS run_demultiplex {} $FASTQDIR $OUTPUT ::: ${bcodefiles[@]}

