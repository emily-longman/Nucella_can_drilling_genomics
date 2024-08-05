#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=genotyping

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --ntasks-per-node=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Call the samtools package
spack load samtools@1.10

#--------------------------------------------------------------------------------

# Working folder
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/mapped_reads

#--------------------------------------------------------------------------------



#--------------------------------------------------------------------------------

# Build bam files
for file in *.sam
    do
    samtools view -t 5 -b $WORKING_FOLDER/$file > $WORKING_FOLDER/${file}.bam
done

# Merge multiple alignment files
samtools merge merged.bam *.bam

# Sort merge bams

# Remove duplicates of final file