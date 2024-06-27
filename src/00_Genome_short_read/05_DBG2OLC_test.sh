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

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# DBG2OLC executable
DBG2OLC=/gpfs1/home/e/l/elongman/software/DBG2OLC

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

# Long reads
#ONT_FOLDER=/netfiles/pespenilab_share/Nucella/raw/ONT
#If you haven't done it yet, unzip the files 
#gunzip $ONT_FOLDER/FC_all.ONT.nuc.fastq.gz
#If you haven't done it yet, unzip the files 
ONT=/netfiles/pespenilab_share/Nucella/raw/ONT/FC_all.ONT.nuc.fastq

#--------------------------------------------------------------------------------
cd $WORKING_FOLDER/SparseAssembler

# Make directory for DBG2OLC
if [ -d "DBG2OLC" ]
then echo "Working DBG2OLC folder exist"; echo "Let's move on."; date
else echo "Working DBG2OLC folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/SparseAssembler/DBG2OLC; date
fi

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
f $ONT