#!/usr/bin/env bash  
#  
#SBATCH -J GL_plink  
#SBATCH -c 4  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Prune for LD

# Load software  
plink=/gpfs1/home/e/l/elongman/software/plink-1.07-x86_64/plink

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
bam_list=$working_folder/info/bam_filelist.list
regions="-rf $working_folder/sites/regions_maf"

#--------------------------------------------------------------------------------

# Parameters for software
WINDOW=25
STEP=10
R=0.5

#--------------------------------------------------------------------------------

$plink --tped $working_folder/plink/Thermofly_all.tped \
--tfam $working_folder/plink/Thermofly_all.tfam \
--indep-pairwise $WINDOW $STEP $R --allow-extra-chr \
--out $working_folder/plink/Thermofly_all.R2.pruned \
--threads 4 