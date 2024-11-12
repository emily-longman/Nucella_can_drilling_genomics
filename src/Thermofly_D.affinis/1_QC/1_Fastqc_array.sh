#!/usr/bin/env bash  
#  
#SBATCH -J Fastqc  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 4:00:00   
#SBATCH --mem 8G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p bluemoon  
#SBATCH --array=239-278%20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

echo ${SLURM_ARRAY_TASK_ID}

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly_D_affinis/
meta=/netfiles/thermofly/METADATA/Thermofly_metadata.vNov11.2024.tsv

#--------------------------------------------------------------------------------

# Use metadata file to extract sample names and forward and reverse reads
FIL1=$(cat ${meta} | awk -F '\t' '{print $24}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
FIL2=$(cat ${meta} | awk -F '\t' '{print $25}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

SAMP_NAME=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo $FIL1
echo $FIL2
echo $SAMP_NAME

#--------------------------------------------------------------------------------

# Create output folder
cd $working_folder
mkdir fastQC

#--------------------------------------------------------------------------------

# Estimate QC statists 
spack load fastqc@0.11.7

fastqc ${FIL1} \
--outdir $working_folder/fastQC

fastqc ${FIL2} \
--outdir $working_folder/fastQC

#--------------------------------------------------------------------------------

# Inform completion of pipeline
echo "done"
date