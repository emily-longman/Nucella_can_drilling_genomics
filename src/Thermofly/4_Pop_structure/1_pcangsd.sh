#!/usr/bin/env bash  
#  
#SBATCH -J pcangsd  
#SBATCH -c 5  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p general  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Generate pca

# Load software  
module load python3.11-anaconda/2024.02-1
module load gcc/13.3.0-xp3epyt
# pcangsd 
export project_dir=/gpfs1/home/e/l/elongman/venv
export venv_name=pcangsd-venv
cd $project_dir
source $venv_name/bin/activate

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
ref=/netfiles/thermofly/GENOMES/basisetae/GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
#bam_list=$working_folder/info/bam_filelist.list
bam_list=$working_folder/info/bam_filelist_reduced.list

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir pcangsd

#--------------------------------------------------------------------------------

# Copy bam list to pcangsd directory
cp $bam_list $working_folder/pcangsd

# Generate covariance matrix with pcangsd
pcangsd \
-b $working_folder/genotype_likelihoods/Thermofly_GL_reduced_minInd_17_depth_6_minMaf_0.1.beagle.gz \
-o $working_folder/pcangsd/Thermofly_SNPs_reduced_minInd_17_depth_6_minMaf_0.1 \
-t 5 

#--------------------------------------------------------------------------------

deactivate