#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=DBG2OLC

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Submit job array
#SBATCH --array=1-16

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------


# SparseAssembler executable
DBG2OLC=/gpfs1/home/e/l/elongman/software/DBG2OLC

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

#--------------------------------------------------------------------------------

# Make directory for DBG2OLC
mkdir $WORKING_FOLDER/SparseAssembler/DBG2OLC

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/SparseAssembler/DBG2OLC

# Use DBG2OLC to construct short but accurate contigs  
$DBG2OLC \
k 51 \
AdaptiveTh 0.0001 \
KmerCovTh 2 \
MinOverlap 20 \
RemoveChimera 1 \
Contigs $WORKING_FOLDER/SparseAssembler/SparseAssembler_test/Contigs.txt \
f ../Pacbio_data/Pacbio_20x.fasta


./DBG2OLC k 17 AdaptiveTh 0.0001 KmerCovTh 2 MinOverlap 20 RemoveChimera 1 Contigs Contigs.txt f ../Pacbio_data/Pacbio_20x.fasta