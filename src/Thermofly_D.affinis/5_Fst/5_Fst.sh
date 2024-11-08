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

# Calculate Fst between groups

# Load software  
spack load angsd@0.933

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6
# Window size for sliding window FST
WINDOW=25000 
# Window step
WINDOW_STEP=5000 

#--------------------------------------------------------------------------------

# Create array
echo ${SLURM_ARRAY_TASK_ID}

array=({1..100})
i="${array[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Calculate Fst for each pairwise comparison

echo " Prepare the fst for easy window analysis etc"

realSFS fst index $working_folder/SFS_sites/Tom/Thermofly_Tom_reduced.saf.idx \
$working_folder/SFS_sites/Olaa/Thermofly_Olaa_reduced_subset_${i}.saf.idx \
-sfs $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites.2dsfs \
-P $CPU -fold 1 -fstout $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites_nMAF


echo "Print SFS priori for each position"
realSFS fst print $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites_nMAF.fst.idx \
-P $CPU > $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites_nMAF_nMAF.bypos.sfs


echo "Get the global estimate of FST throughout the genome"
realSFS fst stats $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites_nMAF.fst.idx \
-P $CPU > $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites_nMAF.fst


echo "Calculate FST by slidingwindow, window size=$WINDOW and step=$WINDOW_STEP"
realSFS  fst stats2 $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites_nMAF.fst.idx \
-win $WINDOW -step $WINDOW_STEP -P $CPU > $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites_nMAF.slidingwindow