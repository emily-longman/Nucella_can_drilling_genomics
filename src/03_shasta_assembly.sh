#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=shasta 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=200G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/shasta.%A_%a.out # Standard output

# Move to the directory where the output files will be saved
cd /netfiles02/pespenilab_share/Nucella/processed/Base_Genome/

#executable
shasta=/gpfs1/home/e/l/elongman/software/shasta/shasta-Linux-0.10.0

#input
infa=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Nuc.3500.fltlong.fastq

#run shasta
$shasta --input $infa --config Nanopore-May2022