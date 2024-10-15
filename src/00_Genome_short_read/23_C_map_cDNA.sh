#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=map_CDNA

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=6-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Request CPU
#SBATCH --cpus-per-task=6

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------


# Index genome 

# Fastqc

# 


STAR=/netfiles/nunezlab/Shared_Resources/Software/STAR-2.7.11b/bin/Linux_x86_64_static





STAR=/netfiles/nunezlab/Shared_Resources/Software/STAR-2.7.11b/bin/Linux_x86_64_static/STAR

WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly
REFERENCE=$WORKING_FOLDER/rename_scaffolds/N.canaliculata_assembly.fasta


# Index genome using STAR
$STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir Nucella_noGFF \
--genomeFastaFiles $REFERENCE 