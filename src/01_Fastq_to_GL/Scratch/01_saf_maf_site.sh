#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=saf_maf_site

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=02-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Request CPU
#SBATCH --cpus-per-task=10

# Submit job array
#SBATCH --array=0-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate saf and maf of SNPs (allele frequencies) for each collection site.
# THIS SCRIPT MATCHES ERICA's and Claire Merot's pipeline. They do it before doing fst, but that makes no sense becuase they don't actually use this input for fst 
# Claire's pipeline says this is used for doing selection analyses later

#--------------------------------------------------------------------------------

#Load modules 
spack load angsd@0.933
spack load samtools@1.10

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
# This is a file with the names . 
arr=("FB" "HC" "MP")
i="${arr[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

# Extract parameters from config file
N_IND=$(wc -l $WORKING_FOLDER/guide_files/${i}_bam.list | cut -d " " -f 1) 
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods_by_site" ]
then echo "Working genotype_likelihoods_by_site folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods_by_site folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods_by_site; date
fi

# Change directory
cd $WORKING_FOLDER/genotype_likelihoods_by_site

if [ -d "${i}" ]
then echo "Working ${i} folder exist"; echo "Let's move on."; date
else echo "Working ${i} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods_by_site/${i}; date
fi

#--------------------------------------------------------------------------------

# Start pipeline

# Specifically calculate the SAF, MAF and GL for all individuals listed in each bam list.

# Move back to working directory
cd $WORKING_FOLDER

# Generate GL's for polymorphic sites for each Nucella collection location
# Since we will use the output for SFS and calculating FSTs/thetas, then we don't want min MAF, p-value filters
angsd \
-b $WORKING_FOLDER/guide_files/${i}_bam.list \
-ref ${REFERENCE} -anc ${REFERENCE} \
-P $NB_CPU \
-nQueueSize 50 \
-doMaf 1 -doSaf 1 -GL 2 -doMajorMinor 3 -doCounts 1 \
-remove_bads 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 \
-minInd $MIN_IND -setMinDepthInd $MIN_DEPTH \
-sites $WORKING_FOLDER/sites_info/sites_all_maf_pruned \
-rf $WORKING_FOLDER/sites_info/regions_all_maf_pruned \
-out $WORKING_FOLDER/genotype_likelihoods_by_site/${i}/${i}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR" 

# -P: number of threads

# -doMaf 1: estimate allele frequencies
# -doSaf 1: estimate the SFS and/or neutrality tests genotype calling
# -GL 2: estimate genotype likelihoods (GL) using the GATK formula
# -doMajorMinor 3: use major and minor from a file
# -doCounts 1: calculate various counts statistics

# -remove_bads 1: remove reads flagged as ‘bad’ by samtools
# -skipTriallelic 1: don’t use sites with >2 alleles
# -uniqueOnly 1: Remove reads that have multiple best hits
# -only_proper_pairs 1: Include only proper pairs (pairs of read with both mates mapped correctly)
# -minMapQ 30: threshold for minimum read mapping quality (Phred)
# -minQ 20: threshold for minimum base quality (Phred)

# -minInd: min number of individuals to keep a site