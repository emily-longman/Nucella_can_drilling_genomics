#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=map_reads

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

# Call the bwa package
module load bwa-0.7.17-gcc-7.3.0-terdbma

#--------------------------------------------------------------------------------

#Define important file locations

#READS indicates the folder where the  reads are stored.
READS=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/clean_Shortreads

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/ShastaRun


#--------------------------------------------------------------------------------


# This indexing step only needs to be done once for the reference file.
bwa index -p ref -a is $REFERENCE/Assembly.fasta

