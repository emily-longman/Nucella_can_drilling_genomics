#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=trim_reads 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script loops through a set of files defined by MYSAMP, matching left and right reads
# and cleans the raw data using fastp according to parameters set below

# cd to the location (path) to the fastq data:
cd /netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads/

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
--html /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/output/fastp/${NAME1}.html \
--json /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/output/fastp/${NAME1}.json

done