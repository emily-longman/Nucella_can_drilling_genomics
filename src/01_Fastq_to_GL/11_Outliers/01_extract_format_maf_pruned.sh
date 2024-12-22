#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=extract_format_maf_pruned

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss
#SBATCH --time=3:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

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

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

# Extract parameters from config file
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1) 
PERC_IND=0.25 # Lower percent ind to 25% for subsequent analyses
MIN_IND_FLOAT=$(echo "($N_IND * $PERC_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "outliers" ]
then echo "Working outliers folder exist"; echo "Let's move on."; date
else echo "Working outliers folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/outliers; date
fi

# Change directory to outliers
cd outliers

if [ -d "baypass" ]
then echo "Working baypass folder exist"; echo "Let's move on."; date
else echo "Working baypass folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/baypass; date
fi

if [ -d "lfmm" ]
then echo "Working lfmm folder exist"; echo "Let's move on."; date
else echo "Working lfmm folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/lfmm; date
fi

if [ -d "rda" ]
then echo "Working rda folder exist"; echo "Let's move on."; date
else echo "Working rda folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/rda; date
fi

#--------------------------------------------------------------------------------

# Change working directory
cd $WORKING_FOLDER

# Input files - all from config 

Rscript $SCRIPT_FOLDER/11_Outliers/01_extract_format_maf_pruned.R "$MIN_MAF" "$PERC_IND" "$MAX_DEPTH_FACTOR"