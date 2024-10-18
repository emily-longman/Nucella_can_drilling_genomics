#!/usr/bin/env bash  
#  
#SBATCH -J Admix_part2  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate admixture part 2

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly

#--------------------------------------------------------------------------------

# Change directory to output
cd $working_folder/ngs_admix

# Loop through the log files to create one output file - identify which value is the lowest for each K - plot that "fopt.gz" for each K
(for log in `ls K_output/*.log`; do grep -Po 'like=\K[^ ]+' $log; done) > $working_folder/ngs_admix/logfile

# K=1 is best