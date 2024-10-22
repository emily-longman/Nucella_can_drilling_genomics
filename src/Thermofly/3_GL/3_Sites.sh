#!/usr/bin/env bash  
#  
#SBATCH -J Sites  
#SBATCH -c 1  
#SBATCH -N 1   
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Get sites list

# Load software  
spack load angsd@0.933
module load R/4.4.0

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir sites

#--------------------------------------------------------------------------------

# Unzip maf file 
#gunzip $working_folder/genotype_likelihoods/Thermofly_GL_reduced_minInd_16_depth_4.mafs.gz

# Change script permissions 
chmod 777 $script_folder/3_GL/3_Sites_list.R

#--------------------------------------------------------------------------------

infile=$working_folder/genotype_likelihoods/Thermofly_GL_reduced_minInd_17_depth_6_minMaf_0.1.mafs
outfile_sites=$working_folder/sites/sites_maf
outfile_regions=$working_folder/sites/regions_maf


Rscript $script_folder/3_GL/3_Sites_list.R "$infile" "$outfile_sites" "$outfile_regions"

angsd sites index $outfile_sites