#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Identify_SNPs

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Submit job array
#SBATCH --array=1-576%20

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

REFERENCE_FAI=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa.fai


#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name of all the bam files to be processed.

SAMPLE_FILE=$WORKING_FOLDER/ALL_bams.txt

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods" ]
then echo "Working genotype_likelihoods folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods; date
fi

#--------------------------------------------------------------------------------

# Calculate genotype likelihoods
angsd -b $SAMPLE_FILE -ref $REFERENCE -out $WORKING_FOLDER/genotype_likelihoods/Nucella \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
-GL 2 -doGlf 4 \
-nThreads 10
# Filter to retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, wihtout trimming, and adjusting for inde/mapping
# -C 50 reduces the effect of reads with excessive mismatches, while -baq 1 computes base alignment quality
# -GL 2: genotype likelihood model as in GATK; -doGlf 4: output in text format


# Genotype calling
angsd -glf $WORKING_FOLDER/genotype_likelihoods/Nucella.glf.gz -fai $REFERENCE_FAI -nInd 15 -out $WORKING_FOLDER/genotype_likelihoods/Nucella \
-doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1

