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

#Load modules 
spack load angsd@0.933

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Path to the directory with the lane merged bams (filtered, sorted and duplicates removed). 
BAMS_FOLDER=$WORKING_FOLDER/bams_merged

#Path to directory with scripts for pipeline
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=40 #change accordingly in SLURM header
echo "using #CPUs ==" $NB_CPU

source $SCRIPT_FOLDER/03_Call_SNPs/00_config.sh

# Unsure what the regions is doing in the angsd script
REGIONS=$WORKING_FOLDER/sites_info/regions_all_maf

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

# Calculate the MAF and GL, with Plink output for LD pruning

# create solo script for bam list

angsd -b $WORKING_FOLDER/genotype_likelihoods_all/Nucella_bam.list \
-ref ${REFERENCE} -anc ${REFERENCE} \
-P $NB_CPU \
-nQueueSize 50 \
-GL 2 -doMajorMinor 1 -doGeno -4 -doPost 1 -postCutoff 0.8 \
-doPlink 2 -doMaf 1 -doCounts 1 \
-remove_bads 1 -skipTriallelic 1 -only_proper_pairs 1 -uniqueOnly 1 /
-minMapQ 30 -minQ 20 \
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 600  \
-out $WORKING_FOLDER/plink/Nucella_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" 

# nQueueSize -50  Maximum number of queud elements