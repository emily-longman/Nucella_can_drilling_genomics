#!/usr/bin/env bash  
#  
#SBATCH -J pcangsd_vis 
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 0:30:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Generate pca using R
module load R/4.4.0

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv
bam_list=$working_folder/info/bam_filelist.list
cov_mat=$working_folder/pcangsd/Thermofly_SNPs.cov

#--------------------------------------------------------------------------------
# Create output folders
cd $working_folder
mkdir figures

#--------------------------------------------------------------------------------

Rscript $script_folder/4_Pop_structure/1_B_pcangsd_vis.R "$cov_mat" "$bam_list" "$meta"


