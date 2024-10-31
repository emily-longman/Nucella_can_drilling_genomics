#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=index_pruned_sites

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=1:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will index the LD-pruned SNPs with angsd.

#--------------------------------------------------------------------------------

#Load modules 
spack load angsd@0.933

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#--------------------------------------------------------------------------------

# Use angsd to index pruned site list

echo "Index sites with angsd"

angsd sites index $WORKING_FOLDER/sites_info/sites_all_maf_pruned