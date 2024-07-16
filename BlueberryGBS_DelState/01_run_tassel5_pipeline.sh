#!/bin/bash

# SLURM parameters
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=08:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 8 processor core(s) per node X 2 threads per core
#SBATCH --mem=48G   # maximum memory per node
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="blueberry_tassel5_discovery_production"
#SBATCH --mail-user=jeffrey.neyhart@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


## 
## TASSEL5 pipeline for Blueberry GBS
## 
## Whole discovery and production pipeline
## 
## This version uses the keyfile corresponding to all unique entries, which have been grouped based
## on ONLY information in the keyfile (i.e. explicitly duplicated sample names.)
## 
## 


# Set error handling options
set -e
set -u
set -o pipefail

# Load the modules
module load tassel5
module load bwa

## Set variables

# Working directory
WD=/project/gifvl_vaccinium/blueberryGenotyping/blueberryEAA

# Name of input directory
INPUT=$WD/input
# Name of input directory with FASTQ
FASTQDIR=$INPUT/fastq_files
# Name of database
DBNAME=$WD/database/blueberry_gbs_discovery.db
# Name of keyfile
KEY=$INPUT/blueberry_gbs_delstate_all_keys.txt

# Name of tag fasta
TAGFASTA=$WD/tags/gbs_tags_for_alignment.fa.gz
# Name of unzipped tag fasta
TAGFASTAUZ=$WD/tags/gbs_tags_for_alignment.fa

# Name of output sam file
SAMOUT=$WD/alignment/gbs_tags_aligned_DraperV1.sam
# Basename of reference index
REFIND=/project/gifvl_vaccinium/bluberryGenotyping/genome_assemblies/V_corymbosum_genome_Draper_v1.0-scaffolds.fasta

# FASTA reference genome
REF=/project/gifvl_vaccinium/bluberryGenotyping/genome_assemblies/V_corymbosum_genome_Draper_v1.0-scaffolds.fasta

# Output stat file
STATSOUTFILE=$WD/stats/snpStats.txt

# VCF Output file
OUTFILE=$WD/snps/blueberryGBS_production_snps.vcf




# Change working directory
cd $WD

## GBSSeqToTagDBPlugin
# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -GBSSeqToTagDBPlugin \
-c 10 \
-db $DBNAME \
-i $INPUT \
-k $KEY \
-e EcoT22I \
-kmerLength 64 \
-minKmerL 20 \
-mnQS 20 \
-mxKmerNum 100000000 \
-endPlugin -runfork1

## TagExportToFastqPlugin 
# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -TagExportToFastqPlugin \
-c 1 \
-db $DBNAME \
-o $TAGFASTA \
-endPlugin -runfork1

# Unzip the fasta file
gunzip -f $TAGFASTA


# Run the alignment
bowtie2 -p 8 --sensitive -x $REFIND -U $TAGFASTAUZ -S $SAMOUT

# # Alternative
# bowtie2 -p 8 --sensitive -x $REFIND1 -U $TAGFASTAUZ -S $SAMOUT1



## SAMToGBSdbPlugin

# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -SAMToGBSdbPlugin  \
-i $SAMOUT \
-db $DBNAME \
-aProp 0.0 \
-aLen 0 \
-minMAPQ 20 \
-endPlugin -runfork1


## DiscoverySNPCallerPluginV2

# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -DiscoverySNPCallerPluginV2 \
-db $DBNAME \
-maxTagsCutSite 1000 \
-mnLCov 0.1 \
-mnMAF 0.00001 \
-ref $REF \
-deleteOldData true \
-endPlugin -runfork1


## SNPQualityProfilerPlugin

# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -SNPQualityProfilerPlugin \
-db $DBNAME \
-statFile $STATSOUTFILE \
-deleteOldData true \
-endPlugin -runfork1


## This took 70 minutes for 13 flowcell-lanes


## 
## ProductionSNPCallerPluginV2
## 
## 


# Execute the plugin
run_pipeline.pl -Xms1G -Xmx48G -fork1 -ProductionSNPCallerPluginV2 \
-db $DBNAME \
-e EcoT22I \
-i $FASTQDIR \
-k $KEY \
-kmerLength 64 \
-do true \
-o $OUTFILE \
-endPlugin -runfork1






