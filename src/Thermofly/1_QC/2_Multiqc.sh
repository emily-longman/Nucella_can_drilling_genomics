#!/usr/bin/env bash  
#  
#SBATCH -J MultiQC 
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 4:00:00   
#SBATCH --mem 8G   
#SBATCH --output=./slurmOutput/%x_%j.out
#SBATCH -p general  

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv

#--------------------------------------------------------------------------------

# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name multiqc #If you haven't already done so, create and name the environment
conda activate multiqc #activate the environment
#conda install -c bioconda multiqc # If you haven't already done so, install the program

#--------------------------------------------------------------------------------

# Create output folder
cd $working_folder
mkdir multiQC

#--------------------------------------------------------------------------------

# Run multiqc on all of the reads
multiqc $working_folder/fastQC \
-n multiqc_report_all.html \
-o multiQC

#--------------------------------------------------------------------------------

echo "done"
date

conda deactivate