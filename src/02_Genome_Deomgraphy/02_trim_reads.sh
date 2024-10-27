#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Map_all_reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Request CPU
#SBATCH --cpus-per-task=10

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Trim reads. 

# Load modules  
fastp=/gpfs1/home/e/l/elongman/software/fastp

#--------------------------------------------------------------------------------

#Define important file locations

# RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/genome_demography

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "trimmed_reads" ]
then echo "Working trimmed_reads folder exist"; echo "Let's move on."; date
else echo "Working trimmed_reads folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/trimmed_reads; date
fi

#--------------------------------------------------------------------------------
