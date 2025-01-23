#!/usr/bin/env bash  
#  
#SBATCH -J Clean_bams  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-22
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Clean bams

# Load software  
module load openjdk/1.8.0
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

echo ${SLURM_ARRAY_TASK_ID}

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.vNov11.2024.tsv

#--------------------------------------------------------------------------------

# Use metadata file to extract sample names
SAMP_NAME=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo $SAMP_NAME

#--------------------------------------------------------------------------------

# Parameters for software
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir bams_qualimap

#--------------------------------------------------------------------------------

# Do QC on cleaned bams
$qualimap bamqc \
-bam $working_folder/bams_clean/${SAMP_NAME}.srt.rmdp.bam \
-outdir $working_folder/bams_qualimap/Qualimap_${SAMP_NAME} \
--java-mem-size=$JAVAMEM