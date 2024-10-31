#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=get_pruned_sites

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=2:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will get a list of the scaffold and positions of LD-pruned SNPs.

#--------------------------------------------------------------------------------

#Load modules 
spack load angsd@0.933
module load R/4.4.0

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file  
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "sites_info" ]
then echo "Working sites_info folder exist"; echo "Let's move on."; date
else echo "Working sites_info folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/sites_info; date
fi

#--------------------------------------------------------------------------------

# Change permissions for script 
chmod 777 $SCRIPT_FOLDER/05_prune_SNPs_in_LD/03_make_sites_list_pruned.R

#--------------------------------------------------------------------------------

# Get list of pruned sites

# Change directory
cd $WORKING_FOLDER/sites_info

# Use R script to format the list of LD-pruned sites

INPUT_plink=$WORKING_FOLDER/plink/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6.R2.pruned.prune.in
INPUT_angsd=$WORKING_FOLDER/sites_info/sites_all_maf

Rscript $SCRIPT_FOLDER/05_prune_SNPs_in_LD/03_make_sites_list_pruned.R "$INPUT_plink" "$INPUT_angsd"

#--------------------------------------------------------------------------------

# Use angsd to index pruned site list

angsd sites index "$INPUT_angsd"_R2.pruned