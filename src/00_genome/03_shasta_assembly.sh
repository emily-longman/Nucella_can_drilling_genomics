#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=shasta 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=200G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/shasta.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

# Move to the directory where the output files will be saved
cd /netfiles/pespenilab_share/Nucella/processed/Base_Genome/

#executable
#shasta=/gpfs1/home/e/l/elongman/software/shasta/shasta-Linux-0.11.1 #get error "No such file or directory", but I followed the installation exactly
shasta=/netfiles/nunezlab/Shared_Resources/Software/shasta/shasta-Linux-0.10.0

# If you haven't done it yet, gunzip the files 
#gunzip /netfiles/pespenilab_share/Nucella/processed/ONT_fltlong/*fastq.gz

#input
infa=/netfiles/pespenilab_share/Nucella/processed/ONT_fltlong/Nuc.3500.fltlong.fastq

#run shasta
$shasta --input $infa --config Nanopore-May2022