#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=plink_LD

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=11

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=2:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

#Load software
plink=/gpfs1/home/e/l/elongman/software/plink-1.07-x86_64/plink

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Define parameters

source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

N_IND=$(wc -l $WORKING_FOLDER/info/Nucella_bam.list | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

WINDOW=25
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

# NOTE: NOT UPDATED YET!
plink --tped angsd/plink/all_plink_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".tped \
--tfam angsd/plink/all_plink_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".tfam \
--indep-pairwise $WINDOW $STEP $R --allow-extra-chr \
--out angsd/plink/all_plink_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".R2.pruned --threads 1 