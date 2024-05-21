#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Genotype_likelihoods

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------


#Load modules 
module load angsd-0.933-gcc-7.3.0-4wsdzjw
module load samtools-1.10-gcc-7.3.0-pdbkohx

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#Folder for joint bams
BAMS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/Merged_Bams

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#--------------------------------------------------------------------------------
# Define parameters

# Java parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods" ]
then echo "Working genotype_likelihoods folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/genotype_likelihoods

#--------------------------------------------------------------------------------

## PREPARE bamlist
# This is a file with the name and full path of all the bam files to be processed.

cd $BAMS_FOLDER
ls -d "$PWD/"* > $OUTPUT/Nucella_bam.list 

#--------------------------------------------------------------------------------

# Estimating Genotype Likelihoods's and allele frequencies for all sites with ANGSD

# File suffix to distinguish analysis choices
SUFFIX_1="GL"

# Generate GL's
angsd -b ${OUTPUT}/Nucella_bam.list \
-ref ${REFERENCE} \
-anc ${REFERENCE} \
-out ${OUTPUT}/Nucella_${SUFFIX_1} \
-nThreads $CPU \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-GL 1 \
-doSaf 1

##### below filters require `do-Counts`
#-doCounts 1 \
#-minInd 4 \
#-setMinDepthInd 1 \
#-setMaxDepthInd 40 \
#-setMinDepth 10 \
#-skipTriallelic 1 \
#-doMajorMinor 1 \
##### below filters require `doMaf`
#-doMaf 1 \
#-SNP_pval 1e-6 \
#-minMaf 0.01

#--------------------------------------------------------------------------------

# Estimate Genotype Likelihoods's and allele frequencies for only the polymorphic sites

# File suffix to distinguish analysis choices
SUFFIX_2="SNPs"

# Generate GL's for polymorphic sites
angsd -b ${OUTPUT}/Nucella_bam.list \
-ref ${REFERENCE} \
-anc ${REFERENCE} \
-out ${OUTPUT}/Nucella_${SUFFIX_2} \
-nThreads $CPU \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-GL 1 \
-doSaf 1 \
-SNP_pval 1e-6