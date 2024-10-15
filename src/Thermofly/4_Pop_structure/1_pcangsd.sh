#!/usr/bin/env bash  
#  
#SBATCH -J pcangsd  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Generate pca

# Load software  
module load python3.11-anaconda
module load gcc/10.5.0
# pcangsd 
export project_dir=/gpfs1/home/e/l/elongman/venv
export venv_name=pcangsd-venv
cd $project_dir
source $venv_name/bin/activate

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked
script_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/Thermofly
bam_list=$working_folder/info/bam_filelist.list

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir pcangsd

#--------------------------------------------------------------------------------

cp $bam_list $working_folder/pcangsd

pcangsd \
-b $working_folder/genotype_likelihoods/Thermofly_GL.beagle.gz \
-o $working_folder/pcangsd/Thermofly_SNPs \
-t 1 

cov_mat=$working_folder/pcangsd/Thermofly_SNPs.cov

#--------------------------------------------------------------------------------

deactivate