#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Get_sites

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will get a list of the scaffold and positions of SNPs identified in the previous script.

# Specifically, it will 
# 1) extract SNP which passed filters
# 2) order the sites file by scaffold names
# 3) make a region file matching the sites files with the same order
# 4) index sites file

#--------------------------------------------------------------------------------

#Load modules 
spack load angsd@0.933
module load R/4.4.0

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#Scripts folder is where all relevant scripts for the pipeline are stored.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Prepare variables 

#Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "sites_info" ]
then echo "Working sites_info folder exist"; echo "Let's move on."; date
else echo "Working sites_info folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/sites_info; date
fi

#--------------------------------------------------------------------------------

# Extract sites
echo "From the maf file, extract a list of SNP chr, position, major all, minor all"

# Unzip maf file if haven't already done so
#gunzip $WORKING_FOLDER/genotype_likelihoods_all/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".mafs.gz

INFILE=$WORKING_FOLDER/genotype_likelihoods_all/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".mafs
OUTFILE_sites=$WORKING_FOLDER/sites_info/sites_all_maf
OUTFILE_regions=$WORKING_FOLDER/sites_info/regions_all_maf

# Change permissions for script 
chmod 777 $SCRIPT_FOLDER/03_Call_SNPs/04_make_sites_list_maxdepth.R

Rscript $SCRIPT_FOLDER/03_Call_SNPs/04_make_sites_list_maxdepth.R "$INFILE" "$OUTFILE_sites" "$OUTFILE_regions"

angsd sites index $OUTFILE_sites