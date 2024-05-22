#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=repeatmasker

# Specify partition
#SBATCH --partition=bluemoon

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

cd /gpfs1/home/e/l/elongman/software

wget https://www.dfam.org/releases/Dfam_3.8/families/FamDB/dfam38_full.0.h5.gz

gunzip dfam38_full.0.h5.gz