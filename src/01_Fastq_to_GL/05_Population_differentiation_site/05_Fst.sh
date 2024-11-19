#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Fst

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate per site Fst (i.e., it will acutally calculate Fst for each pairwise comparison).

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
arr=("FB" "HC" "MP")

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

echo "For sliding window analyses, we will be using a window size of:$WINDOW and a step size of:$WINDOW_STEP"

#--------------------------------------------------------------------------------

# Calculate Fst for the entire genome and using a sliding window

# Number of sites
num_sites="${#arr[@]}" # Length of elements in array

# Estimate Fst for all pairs of collection sites

for i in "${!arr[@]}"; do
for j in "${!arr[@]}"; do
if [ "$i" -lt "$j" ]; then
site1=${arr[i]}
site2=${arr[j]}
echo "Fst between $site1 and $site2"
echo "site 1:" "$site1" 
echo "site 2:" "$site2"

echo "Compute per-site FST indexes"
realSFS fst index \
$WORKING_FOLDER/genotype_likelihoods_by_site/${site1}/${site1}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset.saf.idx \
$WORKING_FOLDER/genotype_likelihoods_by_site/${site2}/${site2}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset.saf.idx \
-sfs $WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset.2dsfs \
-P $NB_CPU -fold 1 -fstout $WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset_nMAF

echo "Print SFS priori for each position"
realSFS fst print \
$WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset_nMAF.fst.idx \
-P $NB_CPU \
> $WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset_nMAF.bypos.sfs
		
echo "Get the global estimate of FST throughout the genome"
realSFS fst stats \
$WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset_nMAF.fst.idx \
-P $NB_CPU \
> $WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset_nMAF.fst
		
echo "Calculate FST by sliding window, window size=$WINDOW and step=$WINDOW_STEP, as given in A_01_config.sh"
realSFS fst stats2 \
$WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset_nMAF.fst.idx \
-win $WINDOW -step $WINDOW_STEP -P $NB_CPU \
> $WORKING_FOLDER/fst/"$site1"_"$site2"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_subset_nMAF.slidingwindow

fi
done
done