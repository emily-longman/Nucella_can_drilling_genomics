#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=gwas_binary_all_SNPs

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=4-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Request CPU
#SBATCH --cpus-per-task=4

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will do a gwas using binary phenotypic data. 
# Based on a generalized linear framework which also allows for quantitative traits and binary and for including additional covariates, using genotype posteriors.

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

# Phenotype data: this file must be a single column with phenotype coded as 0 or 1, each line is one individual in the same order as bamfile.
PHENO=$WORKING_FOLDER/guide_files/Phenotype_Drilled_Binary.txt

# Path to bam list.
BAM_LIST=$WORKING_FOLDER/guide_files/Nucella_bam.list

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=4 #change accordingly in SLURM header
echo "using #CPUs ==" $NB_CPU

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

# Extract parameters from config file
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1) 
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "GWAS" ]
then echo "Working GWAS folder exist"; echo "Let's move on."; date
else echo "Working GWAS folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/GWAS; date
fi

#--------------------------------------------------------------------------------

# Perform GWAS on all samples using binary phenotypic data.

angsd \
-b $BAM_LIST \
-P $NB_CPU \
-nQueueSize 50 \
-yBin $PHENO -doAsso 2 -GL 2 \
-doMaf 1 -doMajorMinor 1 -doCounts 1 -doPost 1 \
-remove_bads 1 -minMapQ 30 -minQ 20 \
-minInd $MIN_IND -setMinDepthInd $MIN_DEPTH -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH \
-sites $WORKING_FOLDER/sites_info/sites_all_maf \
-rf $WORKING_FOLDER/sites_info/regions_all_maf \
-out $WORKING_FOLDER/GWAS/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR".binary.all.SNPs.gwas

# -yBin: File containing binary phenotypic data 
# -doAsso 2: Score Test
# -GL 2: estimate genotype likelihoods (GL) using the GATK formula

# -doMaf 1: estimate allele frequencies
# -doMajorMinor 1: infer the major/minor using different approaches
# -doCounts 1: calculate various counts statistics
# -doPost 1: estimate the posterior genotype probability based on the allele frequency as a prior
# If you use the score statistics -doAsso 2 then calculate the posterior using the allele frequency as prior (-doPost 1).

#--------------------------------------------------------------------------------

# After complete, gunzip the .lrt0.gz output file 