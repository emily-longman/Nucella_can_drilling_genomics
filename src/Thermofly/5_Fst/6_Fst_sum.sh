#!/usr/bin/env bash  
#  
#SBATCH -J Fst_calc  
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

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir Fst_summary

#--------------------------------------------------------------------------------

cd $working_folder/Fst

# Make log of just unweighted 
(for i in 'ls *_allsites_nMAF.fst' ; do cat $i | awk '{print $1}'; done) > $working_folder/Fst_summary/Fst_unweighted_logfile


# Make log of just unweighted 
(for i in 'ls *_allsites_nMAF.fst' ; do cat $i | awk '{print $2}'; done) > $working_folder/Fst_summary/Fst_weighted_logfile