#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=vcf

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=2-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Request CPU
#SBATCH --cpus-per-task=10

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will generate a vcf output. 

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

# Path to bam list.
BAM_LIST=$WORKING_FOLDER/guide_files/Nucella_bam.list

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=10 #change accordingly in SLURM header
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

echo "Calculate the SAF, MAF and GL for all individuals listed in Nucella_bam.list"

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "vcf" ]
then echo "Working vcf folder exist"; echo "Let's move on."; date
else echo "Working vcf folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/vcf; date
fi

#--------------------------------------------------------------------------------

# Generate a vcf output

angsd \
-b $BAM_LIST \
-anc ${REFERENCE} \
-P $NB_CPU \
-nQueueSize 50 \
-doMaf 1 -doSaf 1 -GL 2 -doGlf 2 -doMajorMinor 3 -doGeno 2 -doCounts 1 -doPost 1 -dovcf 1 -doPlink 2 -postCutoff 0.8 \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-minInd 144 -setMaxDepth $MAX_DEPTH \
-sites $WORKING_FOLDER/sites_info/sites_all_maf_pruned \
-rf $WORKING_FOLDER/sites_info/regions_all_maf \
-out $WORKING_FOLDER/vcf/Nucella_SNPs_maf"$MIN_MAF"_pctind"0.75"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pruned

# -doVcf 1: gives a vcf output - simply a wrapper around -gl -dopost -domajorminor -domaf