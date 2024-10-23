#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=bam_list

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=11

# Reserve walltime -- hh:mm:ss
#SBATCH --time=1:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#--------------------------------------------------------------------------------

# Prepare bamlist
# This is a file with the name and full path of all the bam files to be processed.

# Move to bams folder
cd $WORKING_FOLDER/bams_merged

# Create bamlist for all Nucella samples
ls -d "$PWD/"*.bam > $WORKING_FOLDER/guide_files/Nucella_bam.list

#--------------------------------------------------------------------------------

# Make a list of bam files by collection site 

# Create bamlist for HC
ls -d "$PWD/"HC*.bam > $WORKING_FOLDER/guide_files/HC_bam.list

# Create bamlist for FB
ls -d "$PWD/"FB*.bam > $WORKING_FOLDER/guide_files/FB_bam.list

# Create bamlist for MP
ls -d "$PWD/"MP*.bam > $WORKING_FOLDER/guide_files/MP_bam.list
