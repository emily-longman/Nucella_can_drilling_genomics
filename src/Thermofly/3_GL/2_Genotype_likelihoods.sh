#!/usr/bin/env bash  
#  
#SBATCH -J Genotype_likelihoods  
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
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked.fa
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
bam_list=$working_folder/info/bam_filelist_reduced.list
region_file=$working_folder/info/Thermofly_region_file.txt

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
-out $working_folder/genotype_likelihoods/Thermofly_GL_reduced_minInd_14_depth_4_minMaf_0.1 \
-P $CPU \
-doMaf 1 -doSaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-rf $region_file \
-remove_bads 1 -baq 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 \
-minInd 14 -setMinDepthInd 4 -minMaf 0.1 -setMaxDepth 360 \
-SNP_pval 1e-6

#Note: original parameters with all individuals: -minInd 16 -setMinDepthInd 4 -minMaf 0.05 -setMaxDepth 440 \