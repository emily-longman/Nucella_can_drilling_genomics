#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=plink_LD

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Request CPU
#SBATCH --cpus-per-task=5

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Load software
plink=/gpfs1/home/e/l/elongman/software/plink_linux_x86_64_20241022/plink

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

# Path to bam list.
BAM_LIST=$WORKING_FOLDER/guide_files/Nucella_bam.list

#--------------------------------------------------------------------------------

# Define parameters

source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

# Additional parameters for plink
WINDOW_PLINK=25
STEP=10
R=0.5

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "plink" ]
then echo "Working plink folder exist"; echo "Let's move on."; date
else echo "Working plink folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/plink; date
fi

#--------------------------------------------------------------------------------

# Use plink to do LD based pruning 
# (i.e., it will generate pruned subset of markers that are in approximate linkage equilibrium with each other, writing the IDs to plink.prune.in, and the IDs of all excluded variants to plink.prune.out
$plink --tped $WORKING_FOLDER/plink/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6.tped \
--tfam $WORKING_FOLDER/plink/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6.tfam \
--indep-pairwise $WINDOW_PLINK $STEP $R --allow-extra-chr \
--out $WORKING_FOLDER/plink/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6.R2.pruned 

# Variant pruning with indep-pairwise requires 3 parameters: 
# 1) a window size in variant count 
# 2) a variant count to shift the window at the end of each step 
# 3) a pairwise r2 threshold (at each step, pairs of variants in the current window with squared correlation greater than the threshold are noted, and variants are greedily pruned from the window until no such pairs remain)

