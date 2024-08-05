#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Trim_reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=5 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Submit job array
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Trim_reads.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will trim the reads.

# Load modules  
fastp=/gpfs1/home/e/l/elongman/software/fastp

#--------------------------------------------------------------------------------

#Define important file locations

#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#Name of pipeline
PIPELINE=Trim_reads

#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

SAMPLE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/GuideFile_trim_bam.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##               File1                             File2              Snail_ID  Sample#  Lane#    Merged_name    Merged_bam_name
## FB1-1_S84_L002_R1_001.fastq.gz    FB1-1_S84_L002_R2_001.fastq.gz    FB1-1     S84     L002    FB1-1_S84_L002    FB1-1_S84
## FB1-1_S84_L007_R1_001.fastq.gz    FB1-1_S84_L007_R2_001.fastq.gz    FB1-1     S84     L007    FB1-1_S84_L007    FB1-1_S84
## FB1-1_S84_L008_R1_001.fastq.gz    FB1-1_S84_L008_R2_001.fastq.gz    FB1-1     S84     L008    FB1-1_S84_L008    FB1-1_S84
## FB1-2_S173_L002_R1_001.fastq.gz   FB1-2_S173_L002_R2_001.fastq.gz   FB1-2     S173    L002    FB1-2_S173_L002   FB1-2_S173
## ...
## MP9-10_S26_L007_R1_001.fastq.gz   MP9-10_S26_L007_R2_001.fastq.gz   MP9-10    S26     L007    MP9-10_S26_L007   MP9-10_S26
## MP9-10_S26_L008_R1_001.fastq.gz   MP9-10_S26_L008_R2_001.fastq.gz   MP9-10    S26     L008    MP9-10_S26_L008   MP9-10_S26

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $6}' $SAMPLE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
read1=`awk -F "\t" '{print $1}' $SAMPLE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
read2=`awk -F "\t" '{print $2}' $SAMPLE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo $i "+" $read1 "+" $read2

#--------------------------------------------------------------------------------

# This script loops through a set of files defined by MYSAMP, matching left and right reads
# and cleans the raw data using fastp according to parameters set below

# Define the sample code to anlayze

# For each file that has "_R1_001.fastq.gz" (read 1) in the name (the wildcard * allows for the different reps to be captured in the list)
# Start a loop with this file as the input:

for READ1 in *_R1_001.fastq.gz
do

# The partner to this file (read 2) can be found by replacing _R1_001.fastq.gz with _R2_001.fastq.gz
# second part of the input for PE reads

READ2=${READ1/_R1_001.fastq.gz/_R2_001.fastq.gz}

# make the output file names: print the fastq name, replace _# with _#_clean

NAME1=$(echo $READ1 | sed "s/_R1/_R1_clean/g")
NAME2=$(echo $READ2 | sed "s/_R2/_R2_clean/g")

# print the input and output to screen 

echo $READ1 $READ2
echo $NAME1 $NAME2

# call fastp
/gpfs1/home/e/l/elongman/software/fastp -i ${READ1} -I ${READ2} \
-o /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/clean_Shortreads/${NAME1} \
-O /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/clean_Shortreads/${NAME2} \
--detect_adapter_for_pe \
--trim_front1 12 \
--trim_poly_g \
--thread 16 \
--cut_right \
--cut_window_size 6 \
--qualified_quality_phred 20 \
--length_required 35 \
--html /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/processed/fastp/${NAME1}.html \
--json /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/processed/fastp/${NAME1}.json

done