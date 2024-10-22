#!/usr/bin/env bash  
#  
#SBATCH -J SAF_pop_Tom 
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
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked.fa

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir SFS_sites
cd $working_folder/SFS_sites
mkdir Tom

#--------------------------------------------------------------------------------

# Calculate saf for both sites 

#Tom's Trail 
angsd \
-b $working_folder/info/Tom_bam_filelist_reduced.list \
-ref ${ref} -anc ${ref} \
-out $working_folder/SFS_sites/Tom/Thermofly_Tom_reduced \
-P $CPU \
-doMaf 1 -doSaf 1 -GL 2 -doMajorMinor 3 -doCounts 1 \
-sites $working_folder/sites/sites_maf \
-rf $working_folder/sites/regions_maf \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 \
-minInd 6 -setMinDepthInd 6

# Since  we use output for SFS to calculate FSTs/thetas then we don't want min MAF nor p-value filters
# -doMajorMinor 3 means use major and minor from a txt file