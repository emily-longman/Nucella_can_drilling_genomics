#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=partition_genome

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

# Submit job array
#SBATCH --array=1-560

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/partition_genome.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will partition the genome into chunks. 

# It will use the txt file produced from the  previous R script which calculated the number of chunks necessary to split the assembly into 50 contig chunks

# There are 19013 contigs in the assembly. If the assembly is broken into 30 contig chunks then there are 634 chunks. 

#--------------------------------------------------------------------------------
