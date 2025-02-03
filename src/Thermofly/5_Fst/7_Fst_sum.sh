#!/usr/bin/env bash  
#  
#SBATCH -J Fst_sum  
#SBATCH -c 1 
#SBATCH -N 1 # on one node  
#SBATCH -t 1:00:00   
#SBATCH --mem 10G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p general 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Make summary fst document

# Load software  
module load R/4.4.1

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly

#--------------------------------------------------------------------------------

cd $working_folder/Fst_summary


# Call file 
file1=$working_folder/Fst_summary/Fst_unweighted_logfile
file2=$working_folder/Fst_summary/Fst_weighted_logfile

Rscript $script_folder/5_Fst/7_Fst_sum.R "$file1" "$file2"