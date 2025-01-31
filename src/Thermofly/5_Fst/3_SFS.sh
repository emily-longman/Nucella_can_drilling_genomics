#!/usr/bin/env bash  
#  
#SBATCH -J SFS  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out
#SBATCH -p general 
#SBATCH --array=0-99%10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods per group

# Load software  
module load gcc/13.3.0-xp3epyt
module load angsd/0.935-4asngpy

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#--------------------------------------------------------------------------------

# Create array
echo ${SLURM_ARRAY_TASK_ID}

array=({1..100})
i="${array[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir Fst

#--------------------------------------------------------------------------------

# Calculate SFS

realSFS $working_folder/SFS_sites/Tom/Thermofly_Tom_reduced.saf.idx \
$working_folder/SFS_sites/Olaa/Thermofly_Olaa_reduced_subset_${i}.saf.idx \
-P $CPU -maxIter 30 -fold 1 \
> $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites