#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=GL_for_plink

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

# This script will use all bam files to calculate saf, maf and genotype likelihoods on pruned SNP list

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
echo "Keep loci with at least one read for n individuals = $MIN_IND, which is $PERCENT_IND % of total $N_IND individuals"
echo "Filter on allele frequency = $MIN_MAF"

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods_all_pruned" ]
then echo "Working genotype_likelihoods_all_pruned folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods_all_pruned folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods_all_pruned; date
fi

#--------------------------------------------------------------------------------

# Get list of pruned sites

# Use R script to format the list of LD-pruned sites, then have angsd index it

INPUT_plink=$WORKING_FOLDER/plink/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6.R2.pruned.prune.in
INPUT_angsd=$WORKING_FOLDER/sites_info/sites_all_maf

#Rscript Rscripts/make_site_list_pruned.r "$INPUT_plink" "$INPUT_angsd"
angsd sites index "$INPUT_angsd"_R2.pruned

#--------------------------------------------------------------------------------

# Calculate the MAF and GL, with Plink output for LD pruning

# create solo script for bam list

angsd \
-b $BAM_LIST \
-ref ${REFERENCE} -anc ${REFERENCE} \
-P $NB_CPU \
-nQueueSize 50 \
-GL 2 -doMajorMinor 1 -doGeno -4 -doPost 1 -postCutoff 0.8 \
-doPlink 2 -doMaf 1 -doCounts 1 \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 $REGIONS \
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH \
-out $WORKING_FOLDER/plink/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6 

# nQueueSize -50  Maximum number of queud elements
# Notice the extra minus in the -dogeno -4 argument, this will suppress the -doGeno output.
# postCutoff: Call only a genotype with a posterior above this threshold.