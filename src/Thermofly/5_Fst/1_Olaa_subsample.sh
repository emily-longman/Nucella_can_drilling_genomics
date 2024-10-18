#!/usr/bin/env bash  
#  
#SBATCH -J Olaa_array 
#SBATCH -c 1 
#SBATCH -N 1 # on one node  
#SBATCH -t 2:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Randomly sample 7 out of 11 individuals from Olaa Forest

# Load software  
module load R/4.4.0

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder/info
mkdir Olaa_subsample
cd $working_folder/info/Olaa_subsample

#--------------------------------------------------------------------------------

infile=$working_folder/info/Olaa_bam_filelist_reduced.list

Rscript $script_folder/5_Fst/1_Olaa_subsample.R "$infile" 