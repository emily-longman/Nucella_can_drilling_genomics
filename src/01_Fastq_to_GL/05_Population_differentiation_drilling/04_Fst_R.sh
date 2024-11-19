#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Fst_R_drilling

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=20:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=1G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will do the second step in calculating the 2d site frequency spectrums (sfs) and Fst for each collection site pair.

#--------------------------------------------------------------------------------

#Load modules 
module load R/4.4.0

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_mask/N.canaliculata_assembly.fasta.masked

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Establish the array
arr=("Drilled" "Not.Drilled")

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

#--------------------------------------------------------------------------------

# Sum the sfs for each site pair

# Number of sites
num_sites="${#arr[@]}" # Length of elements in array

# Run pairwise FST R script for all populations listed

for i in "${!arr[@]}"; do
for j in "${!arr[@]}"; do
if [ "$i" -lt "$j" ]; then
group1=${arr[i]}
group2=${arr[j]}
echo "Fst between $group1 and $group2"
echo "group 1:" "$group1" 
echo "group 2:" "$group2"

file=$WORKING_FOLDER/fst_drilling/"$group1"_"$group2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset

Rscript $SCRIPT_FOLDER/05_Population_differentiation_drilling/04_Fst.R "$file"

fi
done
done

# Note: Change permissions on .saf.idx files

# -maxIter: maximum number of iterations in the EM algorithm
# -fold 1: estimate the folded SFS (need to do if you don't have an ancestral state)
