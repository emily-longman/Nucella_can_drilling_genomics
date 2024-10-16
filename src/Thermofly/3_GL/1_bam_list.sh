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

# Change working directory
cd $working_folder/bams_clean

# Create bamlist
ls -d "$PWD/"*.bam > $working_folder/info/bam_filelist.list

# Create bamlist for Tom's Trail
ls -d "$PWD/"*B_T*.bam > $working_folder/info/Tom_bam.list

# Create bamlist for Olaa Forest
ls -d "$PWD/"*B_O*.bam > $working_folder/info/Olaa_bam.list

# Change working directory
cd $working_folder/bams_clean_reduced

# Create bamlist (removing samples w/ coverage <4)
ls -d "$PWD/"*.bam > $working_folder/info/bam_filelist_reduced.list