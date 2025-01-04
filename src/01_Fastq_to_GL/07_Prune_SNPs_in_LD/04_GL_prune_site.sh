#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=GL_pruned_site

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

# Submit job array
#SBATCH --array=0-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will use all bam files to calculate saf, maf and genotype likelihoods on pruned SNP list for drilled and non-drilled snails

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
i="${arr[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

# Extract parameters from config file
N_IND=$(wc -l $WORKING_FOLDER/guide_files/${i}_bam.list | cut -d " " -f 1) 
PERC_IND=0.25 # Lower percent ind to 25% for subsequent analyses
MIN_IND_FLOAT=$(echo "($N_IND * $PERC_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods_site_pruned" ]
then echo "Working genotype_likelihoods_site_pruned folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods_site_pruned folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods_site_pruned; date
fi

#--------------------------------------------------------------------------------

# Calculate the MAF and GL, with Plink output for LD pruning for each drilling group

# Move back to working directory
cd $WORKING_FOLDER

echo "Working on group ${i}, with $N_IND individuals."
echo "Will filter for sites with at least one read in $MIN_IND individuals, which is $PERC_IND of the total."

# Generate GL's for polymorphic sites for each Nucella drilling group 
angsd \
-b $WORKING_FOLDER/guide_files/${i}_bam.list \
-ref ${REFERENCE} -anc ${REFERENCE} \
-P $NB_CPU \
-nQueueSize 50 \
-doMaf 1 -doSaf 1 -GL 2 -doGlf 2 -doMajorMinor 3 -doGeno 2 -doPost 1 -doBcf 1 \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 \
-minInd $MIN_IND \
-sites $WORKING_FOLDER/sites_info/sites_all_maf_pruned \
-rf $WORKING_FOLDER/sites_info/regions_all_maf \
-out $WORKING_FOLDER/genotype_likelihoods_site_pruned/${i}_SNPs_maf"$MIN_MAF"_pctind"$PERC_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6_pruned

# Note reduced minInd at this stage to 25%