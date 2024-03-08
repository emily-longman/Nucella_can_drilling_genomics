#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fltlong 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G #<= this may depend on your resources

# Submit job array
#SBATCH --array=0-3

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/fltlong.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to


filtlong=/gpfs1/home/e/l/elongman/software/Filtlong/bin/filtlong
input=/netfiles/pespenilab_share/Nucella/raw/ONT/FC_all.ONT.nuc.fastq.gz

arr=(1000 2000 3500 5000)
L="${arr[$SLURM_ARRAY_TASK_ID]}"
echo $L

$filtlong \
--min_length $L \
$input | gzip > Nuc.$L.fltlong.fastq.gz