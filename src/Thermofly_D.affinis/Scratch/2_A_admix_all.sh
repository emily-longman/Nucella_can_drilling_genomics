#!/usr/bin/env bash  
#  
#SBATCH -J Admix  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Calculate admixture

# Load software  
NGSadmix=/gpfs1/home/e/l/elongman/software/NGSadmix

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
region_file=$working_folder/info/Thermofly_region_file.txt
bam_list=$working_folder/info/bam_filelist_reduced.list

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6

#Min number of pop to consider for NGS admix
K_MIN=1
#Maximum number of pop to consider for NGS admix
K_MAX=4

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir ngs_admix
cd $working_folder/ngs_admix
mkdir K_output

#--------------------------------------------------------------------------------

## Run admixture for Kmin to Kmax 10 times each
for j in {1..10}; 
do

for i in $(seq $K_MIN $K_MAX)
do 
echo $i
$NGSadmix -P $CPU \
-likes $working_folder/genotype_likelihoods/Thermofly_GL_reduced_minInd_16_depth_4.beagle.gz \
-minMaf $0.05 -K $i -o $working_folder/ngs_admix/K_output/Thermofly_GL_reduced_minInd_16_depth_4_maf_K{$i}_run{$j}
done
done