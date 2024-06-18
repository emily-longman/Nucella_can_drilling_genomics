#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Genotype_likelihoods_groups

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Submit job array
#SBATCH --array=0-1

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/GL_sites.%A_%a.out # Standard output

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

#Name of pipeline
PIPELINE=Genotype_Likelihoods_groups

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#--------------------------------------------------------------------------------

# Establish array

arr=("Drilled" "NotDrilled")
L="${arr[$SLURM_ARRAY_TASK_ID]}"
echo $L

#--------------------------------------------------------------------------------

## Bamlists
# These are files with the name and full path of all the bam files to be processed for snails that both Drilled and NotDrilled

#Folder for joint bams
BAM_LIST_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/Guide_Files

#--------------------------------------------------------------------------------


# Move to working directory
cd $WORKING_FOLDER

# Begin Pipeline

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods" ]
then echo "Working genotype_likelihoods folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/genotype_likelihoods

#--------------------------------------------------------------------------------

# Estimating Genotype Likelihoods's and allele frequencies for all sites with ANGSD

# Move back to working directory
cd $WORKING_FOLDER

# File suffix to distinguish analysis choices
SUFFIX_1="GL_allsites"

# Generate GL's for each drilling group
#angsd -b ${BAM_LIST_FOLDER}/${L}_bam.list \
#-ref ${REFERENCE} \
#-anc ${REFERENCE} \
#-out ${OUTPUT}/${L}_${SUFFIX_1} \
#-nThreads $CPU \
#-remove_bads 1 \
#-C 50 \
#-baq 1 \
#-minMapQ 30 \
#-minQ 20 \
#-skipTriallelic 1 \
#-GL 1 \
#-doCounts 1 \
#-doMajorMinor 1 \
#-doMaf 1 \
#-doSaf 1 \
#-doHWE 1


# Generate GL's for each group
angsd -b ${BAM_LIST_FOLDER}/${L}_bam.list \
-ref ${REFERENCE} \
-anc ${REFERENCE} \
-out ${OUTPUT}/${L}_${SUFFIX_1} \
-nThreads $CPU \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-GL 1 \
-doSaf 1
#--------------------------------------------------------------------------------

# Inform that the pipeline is done

echo "pipeline completed" $(date)