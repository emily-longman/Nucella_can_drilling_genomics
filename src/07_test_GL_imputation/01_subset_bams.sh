#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Subset_bams

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
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------


#Load modules 
module load samtools-1.10-gcc-7.3.0-pdbkohx

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#Folder for joint bams
BAMS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/Merged_Bams

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#--------------------------------------------------------------------------------

cd $WORKING_FOLDER

# Make directory for each DBG2OLC parameter combination
if [ -d "BAMS_subset" ]
then echo "Working BAMS_subset folder exist"; echo "Let's move on."; date
else echo "Working BAMS_subset folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/BAMS_subset; date
fi

#--------------------------------------------------------------------------------

# Change to output folder
cd $WORKING_FOLDER/BAMS_subset

# Samtools coverage to determine scaffold coverage

#Make temporary linefile with list of input BAM files
#ls $WORKING_FOLDER/Merged_Bams/*.Lanes_merged.bam > guide.txt

#samtools coverage -m -b guide.txt

#--------------------------------------------------------------------------------

# Change to output folder
cd $WORKING_FOLDER/BAMS_subset

# Examples form online
#samtools view -b input.bam "Chr10:18000-45500" > output.bam



# Testing 
samtools view -b $WORKING_FOLDER/Merged_Bams/${i}.Lanes_merged.bam "ntLink_33941:ntLink_33950" > $WORKING_FOLDER/BAMS_subset/${i}.subset.bam

