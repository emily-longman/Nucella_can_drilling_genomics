#!/usr/bin/env bash  
#  
#SBATCH -J SFS_pop  
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
region_file=$working_folder/info/Thermofly_region_file.txt

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir SFS_sites

#--------------------------------------------------------------------------------

# Calculate SFS for both sites - need different parameters since different number of ind in each

#Tom's Trail (7 ind)
angsd \
-b $working_folder/info/Tom_bam.list \
-ref ${ref} -anc ${ref} \
-out $working_folder/SFS_sites/Thermofly_Tom_saf \
-P $CPU \
-doMaf 1 -doSaf 1 -GL 2 -doMajorMinor 1 \
-rf $region_file \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-minInd 5 -setMinDepthInd 4 -setMaxDepth 140 \

#Olaa Forest (15 ind)
angsd \
-b $working_folder/info/Olaa_bam.list \
-ref ${ref} -anc ${ref} \
-out $working_folder/genotype_likelihoods_groups/Thermofly_Olaa \
-P $CPU \
-doMaf 1 -doSaf 1 -GL 2 -doMajorMinor 3 -doCounts 1 \
-sites $working_folder/sites/sites_maf \
-rf $working_folder/sites/regions_maf \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-minInd 11 -setMinDepthInd 4 -setMaxDepth 300 \

# Since  we use output for SFS to calculate FSTs/thetas then we don't want min MAF nor p-value filters