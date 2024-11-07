#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ngsF_inbred

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Request CPU
#SBATCH --cpus-per-task=4

# Submit job array
#SBATCH --array=0-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate inbreeding coefficients with ngsF (https://github.com/fgvieira/ngsF/tree/master)

#--------------------------------------------------------------------------------

#Load modules 
module load singularity/3.7.1

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_mask/N.canaliculata_assembly.fasta.masked

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Move to working directory
cd $SCRIPT_FOLDER/06_Inbreeding_coefficients

# If you haven't already done so, get sif 
#singularity pull shub://seedpcseed/metaerg:latest

#--------------------------------------------------------------------------------

# Export path to braker sif
export ngsTools_SIF=$SCRIPT_FOLDER/06_Inbreeding_coefficients/metaerg_latest.sif

# Change permissions for associated script
chmod +x 02_ngsF_inbred.sh

# Use the singularity to run angsd
singularity exec $SCRIPT_FOLDER/06_Inbreeding_coefficients/metaerg_latest.sif $SCRIPT_FOLDER/06_Inbreeding_coefficients/02_ngsF_inbred.sh