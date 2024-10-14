#!/usr/bin/env bash  
#  
#SBATCH -J bam_list  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 1:00:00   
#SBATCH --mem 10G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir info

#--------------------------------------------------------------------------------

# Create bamlist
cd $working_folder/bams_clean
ls -d "$PWD/"*.bam > $working_folder/info/bam_filelist.list