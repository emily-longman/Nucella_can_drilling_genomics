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


# I don't get what is trying to be extracted here...
# I get we are looping through the log files and searching using grep for a pattern using perl, but what pattern are we looking for? Maybe my output is incorrect and that is the problem?
(for log in `ls K_output/*.log`; do grep -Po 'like=\K[^ ]+' $log; done) > $working_folder/ngs_admix/logfile