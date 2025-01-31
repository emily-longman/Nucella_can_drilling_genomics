#!/usr/bin/env bash  
#  
#SBATCH -J SAF_pop_Olaa_array 
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p general 
#SBATCH --array=0-99%10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods per group

# Load software  
module load gcc/13.3.0-xp3epyt
module load angsd/0.935-4asngpy

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.vNov11.2024.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#--------------------------------------------------------------------------------

# Create array
echo ${SLURM_ARRAY_TASK_ID}

array=({1..100})
i="${array[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir SFS_sites
cd $working_folder/SFS_sites
mkdir Olaa

#--------------------------------------------------------------------------------

# Calculate saf for both sites 

#Olaa Forest
angsd \
-b $working_folder/info/Olaa_subsample/Olaa_bam_filelist_reduced_${i}.list \
-ref ${ref} -anc ${ref} \
-out $working_folder/SFS_sites/Olaa/Thermofly_Olaa_reduced_subset_${i} \
-P $CPU \
-doMaf 1 -doSaf 1 -GL 2 -doMajorMinor 1 \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 \
-minInd 6 -setMinDepthInd 6


# Since  we use output for SFS to calculate FSTs/thetas then we don't want min MAF nor p-value filters
# -doMajorMinor 3 means use major and minor from a txt file
