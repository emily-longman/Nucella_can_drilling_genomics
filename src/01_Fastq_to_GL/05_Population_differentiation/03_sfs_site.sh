#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=sfs_sites

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=02-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G 

# Request CPU
#SBATCH --cpus-per-task=10

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate the 2d site frequency spectrums (sfs) and Fst for each collection site pair.

#--------------------------------------------------------------------------------

#Load modules 
spack load angsd@0.933

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_mask/N.canaliculata_assembly.fasta.masked

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=10 #change accordingly in SLURM header
echo "using #CPUs ==" $NB_CPU

#--------------------------------------------------------------------------------

# Establish the array
# This is a file with the names of the collection sites. 
arr=("FB" "HC" "MP")

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "fst" ]
then echo "Working fst folder exist"; echo "Let's move on."; date
else echo "Working fst folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fst; date
fi

#--------------------------------------------------------------------------------

# Create site frequency spectrums for each site pair

# Number of sites
num_sites="${#arr[@]}" # Length of elements in array

# Estimate pairwise FST for all populations listed

for i in "${!arr[@]}"; do
for j in "${!arr[@]}"; do
if [ "$i" -lt "$j" ]; then
site1=${arr[i]}
site2=${arr[j]}
echo "Fst between $site1 and $site2"
echo "site 1:" "$site1" 
echo "site 2:" "$site2"

echo "Calculate the 2dsfs priors"

realSFS \
$WORKING_FOLDER/genotype_likelihoods_by_site/${site1}/${site1}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset.saf.idx \
$WORKING_FOLDER/genotype_likelihoods_by_site/${site2}/${site2}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset.saf.idx \
-P $NB_CPU -maxIter 30 -fold 1 \
> $WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset

fi
done
done

# Note: Change permissions on .saf.idx files

# -maxIter: maximum number of iterations in the EM algorithm
# -fold 1: estimate the folded SFS (need to do if you don't have an ancestral state)