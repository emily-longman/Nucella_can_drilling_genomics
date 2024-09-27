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

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=1:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#Path to the directory with the lane merged bams (filtered, sorted and duplicates removed). 
BAMS_FOLDER=$WORKING_FOLDER/bams_merged

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "info" ]
then echo "Working info folder exist"; echo "Let's move on."; date
else echo "Working info folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/info; date
fi

#--------------------------------------------------------------------------------

# Prepare bamlist
# This is a file with the name and full path of all the bam files to be processed.

# Move to bams folder
cd $BAMS_FOLDER

# Create bamlist for all Nucella samples
ls -d "$PWD/"* > $WORKING_FOLDER/info/Nucella_bam.list