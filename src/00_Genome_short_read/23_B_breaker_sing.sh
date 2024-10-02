#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=breaker_sing

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=2

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

# Load modules
module load singularity/3.7.1

#Define important file locations

#Working folder is core folder where this pipeline is being run.
SCRIPTS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read

#--------------------------------------------------------------------------------

# Change directory
cd $SCRIPTS_FOLDER

# Export path to braker sif
export BRAKER_SIF=$SCRIPTS_FOLDER/braker3.sif

# Change permissions for associated script
chmod +x 23_breaker.sh

# Use the singularity to run angsd
singularity exec $SCRIPTS_FOLDER/braker3.sif $SCRIPTS_FOLDER/23_breaker.sh