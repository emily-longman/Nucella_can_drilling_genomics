#!/usr/bin/env bash  
#  
#SBATCH -J GL_plink  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods

# Load software  
spack load angsd@0.933

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
CPU=6

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir plink

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods
angsd \
-b $bam_list \
-ref ${ref} -anc ${ref} \
-out $working_folder/plink/Thermofly_all \
-P $CPU \
-nQueueSize 50 \
-doMaf 1 -GL 2 -doMajorMinor 1 -doCounts 1 -doGeno -4 -doPost 1 -postCutoff 0.8 -doPlink 2 \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 $regions \
-minInd 16 -minMaf 0.05 -setMaxDepth 352 \
-SNP_pval 1e-6