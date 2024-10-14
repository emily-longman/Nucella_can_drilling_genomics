#!/usr/bin/env bash  
#  
#SBATCH -J Genotype_likelihoods  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon  
#SBATCH --array=1-22
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Clean bams

# Load software  
spack load angsd@0.933

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
bam_list=$working_folder/info/bam_filelist.list

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir genotype_likelihoods

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods
angsd \
-b $bam_list \
-ref ${ref} -anc ${ref} \
-out $working_folder/genotype_likelihoods/Thermofly_GL \
-P $CPU \
-doMaf 1 -doSaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-remove_bads 1 -baq 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 \
-minInd 16 -setMinDepthInd 1 -minMaf 0.05 -setMaxDepth 352 \
-SNP_pval 1e-6