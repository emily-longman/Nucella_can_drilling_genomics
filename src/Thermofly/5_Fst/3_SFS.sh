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

# Calculate genotype likelihoods per group

# Load software  
spack load angsd@0.933

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir Fst

#--------------------------------------------------------------------------------

# Calculate SFS

realSFS $working_folder/SFS_sites/Tom/Thermofly_Tom_reduced.saf.idx \
$working_folder/SFS_sites/Olaa/Thermofly_Olaa_reduced_subset_${i}.saf.idx \
-P $CPU -maxIter 30 -fold 1 \
> $working_folder/Fst/Thermofly_Tom_Olaa_subset_${i}_allsites