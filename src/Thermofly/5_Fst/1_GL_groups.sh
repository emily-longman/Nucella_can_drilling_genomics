#!/usr/bin/env bash  
#  
#SBATCH -J Genotype_likelihoods  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/GL_groups.%A_%a.out 
#SBATCH -p bluemoon 
# Submit job array
#SBATCH --array=0-1 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods per group

# Load software  
spack load angsd@0.933

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked
#region_file=$working_folder/info/Thermofly_region_file.txt

#--------------------------------------------------------------------------------

#  Array
arr=("Tom" "Olaa")
L="${arr[$SLURM_ARRAY_TASK_ID]}"
echo $L

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir genotype_likelihoods_groups

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods
angsd \
-b $working_folder/info/${L}_bam.list \
-ref ${ref} -anc ${ref} \
-out $working_folder/genotype_likelihoods_groups/Thermofly_${L} \
-P $CPU \
-doMaf 1 -doSaf 1 -GL 2 -doMajorMinor 3 -doCounts 1 \
-sites $working_folder/sites/sites_maf /
-rf $working_folder/sites/regions_maf \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 \
-minInd 16 -setMinDepthInd 4 -setMaxDepth 440 \

# Since  we use output for SFS to calculate FSTs/thetas then we don't want min MAF nor p-value filters