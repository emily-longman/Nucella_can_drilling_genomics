#!/usr/bin/env bash  
#  
#SBATCH -J Fst_R  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon 
#SBATCH --array=0-99%10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate Fst between groups

# Load software  
module load R/4.4.0

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly

#--------------------------------------------------------------------------------


# Create array
echo ${SLURM_ARRAY_TASK_ID}

array=({1..100})
i="${array[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Call file 
file=$working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites

Rscript $script_folder/5_Fst/4_Fst.R "$file"