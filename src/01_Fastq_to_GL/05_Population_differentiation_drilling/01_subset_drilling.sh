#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=subset_drilling

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=2G 

# Submit job array
#SBATCH --array=0-1

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will randomly subset the sites such that each collection site has the same number.
# Differing number of individuals among populations can strongly influence Fst values.

#--------------------------------------------------------------------------------

#Load modules 
module load R/4.4.0

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Establish the array  
arr=("Drilled" "Not.Drilled")
i="${arr[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Change working directory
cd $WORKING_FOLDER/guide_files

# Randomly subset bam file lists
# Drilled has 75 individuals, while Not Drilled has 117 - thus subset the Not Drilled bamlists to also have 75 individuals

# Input files
Drilling_bamlist=$WORKING_FOLDER/guide_files/Nucella_bam_${i}.list
Drilling_group=${i}

Rscript $SCRIPT_FOLDER/05_Population_differentiation_drilling/01_subset_drilling.R "$Drilling_bamlist" "$Drilling_group" 