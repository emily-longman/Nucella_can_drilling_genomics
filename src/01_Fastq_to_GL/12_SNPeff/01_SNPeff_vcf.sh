#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SNPeff_vcf

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run SNPeff on a vcf

#Load modules 
module load snpeff

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# vcf is the vcf input file
vcf=/gpfs1/home/e/l/elongman/scratch/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/vcf_pruned/Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2_pval1e6_pruned.vcf

# datadir is the directory where SNPeff directory was built
datdir=/netfiles/pespenilab_share/Nucella/processed

# param is the SNPeff config file
param=/netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/snpEff.config

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "SNPeff" ]
then echo "Working SNPeff folder exist"; echo "Let's move on."; date
else echo "Working SNPeff folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/SNPeff; date
fi

#--------------------------------------------------------------------------------

# Run SNPeff on the vcf file

# Change directory
cd $WORKING_FOLDER/SNPeff

# Run SNPeff
snpeff -c $param -dataDir $datdir N.can_genome_Dec2024 $vcf > $WORKING_FOLDER/SNPeff/Nucella_SNPs.pruned.annotate.vcf