#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Merge_bams

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will merge all of the bam files together. 

# Load modules  
spack load samtools@1.10
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/genome_demography

# GL folder is folder where the bam.
BAM_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/bams_clean

#-------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "bam" ]
then echo "Working bam folder exist"; echo "Let's move on."; date
else echo "Working bam folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bam; date
fi

if [ -d "bam_qualimap" ]
then echo "Working bam_qualimap folder exist"; echo "Let's move on."; date
else echo "Working bam_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bam_qualimap; date
fi

#--------------------------------------------------------------------------------

# I will merge all of the bam files produced in the fastq to GL pipeline into one bam file

# Make temporary linefile with list of input BAM files
ls $BAM_FOLDER/*.srt.rmdp.bam > $WORKING_FOLDER/N.can.guide.txt

# Merge the 3 sequencing lanes
samtools merge \
-b $WORKING_FOLDER/N.can.guide.txt \
$WORKING_FOLDER/bam/N.can.bam

# Remove the temporary guide file
rm $WORKING_FOLDER/N.can.guide.txt

# Assess quality of final file
$qualimap bamqc \
-bam $WORKING_FOLDER/bam/N.can.bam \
-outdir $WORKING_FOLDER/bam_qualimap \
--java-mem-size=$JAVAMEM

#--------------------------------------------------------------------------------

# Index bams with samtools
samtools index $WORKING_FOLDER/bam/N.can.bam

#--------------------------------------------------------------------------------