#!/usr/bin/env bash  
#  
#SBATCH -J SFS  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon 
# Submit job array
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate Fst between groups

# Load software  
spack load angsd@0.933

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
file=$working_folder/Fst/Thermofly_Tom_Olaa

#--------------------------------------------------------------------------------

realSFS fst index \
$working_folder/SFS_sites/Thermofly_Tom.saf.idx \
$working_folder/SFS_sites/Thermofly_Olaa.saf.idx \
-sfs $working_folder/Fst/Thermofly_Tom_Olaa_allsites.sfs \
-fold 1 -fstout $working_folder/SFS_sites/Tom_Olaa_allsites -whichFst 1