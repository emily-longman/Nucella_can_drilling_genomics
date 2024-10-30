#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=qualimap_multi_bamqc

# Specify partition
#SBATCH --partition=bigmem

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=500G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------
# This script will calculate average coverage 

#Load modules 
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#--------------------------------------------------------------------------------
# Define parameters
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Move to working folder
cd $WORKING_FOLDER

# Generate output folder 
if [ -d "qualimap_multi_bamqc" ]
then
echo "Working qualimap_multi_bamqc folder exist"
echo "lets move on"
else 
echo "Folder doesnt exist. Lets fix that"
mkdir qualimap_multi_bamqc
fi

OUTPUT=$WORKING_FOLDER/qualimap_multi_bamqc

#--------------------------------------------------------------------------------

## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

GUIDE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/Guide_Files/Qualimap_bam_list.txt

# Use script below to make list of the path to the folder which contains results of BAM QC analysis. 
# Download and open file in excel. Add column with name of the sample
# Move to qualimap folder
#cd $WORKING_FOLDER/Merged_Bams_qualimap
# Create a list for each collection location
#ls -d "$PWD/"Qualimap_LaneMerged_* > $OUTPUT/Qualimap_bam_list.txt 

#--------------------------------------------------------------------------------

# Assess quality of all bam files
$qualimap multi-bamqc \
-d $GUIDE_FILE \
-outdir $OUTPUT \
--java-mem-size=$JAVAMEM