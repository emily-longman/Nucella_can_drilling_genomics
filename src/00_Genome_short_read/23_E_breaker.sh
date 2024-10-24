#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=breaker

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=4

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location where the reference genome. (note: copied a final version from pilon to repeatmasker directory)
REFERENCE=$WORKING_FOLDER_SCRATCH/repeatmasker/polished_assembly.fasta
#Will need to update...

#Working folder is core folder where this pipeline is being run.
SCRIPTS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read

# Export path to braker sif
export BRAKER_SIF=$SCRIPTS_FOLDER/braker3.sif

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "breaker" ]
then echo "Working breaker folder exist"; echo "Let's move on."; date
else echo "Working breaker folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/breaker; date
fi

#--------------------------------------------------------------------------------

# Execute breaker
braker.pl --species=Nucella_canaliculata --genome=$REFERENCE \
--bam=file1.bam
