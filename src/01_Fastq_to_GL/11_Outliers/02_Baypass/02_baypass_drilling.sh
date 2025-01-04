#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_2_drilling

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10

# Reserve walltime -- hh:mm:ss
#SBATCH --time=3:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run baypass on the snail data. 
# It will use the gfile produced in step 01_reformat and the omega file produced in step 02_Baypass/01_baypass.sh  

#Load modules 
module load gcc/10.5.0
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

# Path to bam list.
BAM_LIST=$WORKING_FOLDER/guide_files/Nucella_bam.list

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

# Extract parameters from config file
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1) 
PERC_IND=0.25 # Lower percent ind to 25% for subsequent analyses
MIN_IND_FLOAT=$(echo "($N_IND * $PERC_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

# Number of groups/populations
npop=2

#--------------------------------------------------------------------------------

# Move to baypass directory
cd $WORKING_FOLDER/outliers/baypass

# Run baypass using omega file 
$baypass -npop $npop \
-gfile $WORKING_FOLDER/outliers/baypass/by_group_"$MIN_MAF"_pctind"$PERC_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR".mafs.pruned.baypass \
-omegafile $WORKING_FOLDER/outliers/baypass/Drilling_r2_mat_omega.out \
-outprefix Drilling_r2.output -nthreads 10