#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=breaker_array

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=1-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G

# Submit job array
#SBATCH --array=1 #1-631%30

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will annotate the soft masked genome. However, given the large size of the genome, the length of time exceeds run limits on the VACC. 
# Thus, to make the process more manageable, the genome is broken up into individual scaffolds then braker is run on each scaffold.
# To do this utilize both an array and a while loop. 
# The previous script produced a guide file that groups scaffolds into 30 scaffold chunks, for a total of 631 partitions.
# For each partition, this script will loop over each scaffold name, break the genome and bam file into that scaffold then clean that scaffold.

# Load modules
module load singularity/3.7.1

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location of the soft-masked reference genome.
REFERENCE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/N.canaliculata_assembly.fasta.softmasked

#Working folder is core folder where this pipeline is being run.
SCRIPTS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read

# Export path to braker sif
export BRAKER_SIF=$SCRIPTS_FOLDER/23_braker_singularity/braker3.sif

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/braker_ab_initio

if [ -d "braker_ab_initio" ]
then echo "Working braker_ab_initio folder exist"; echo "Let's move on."; date
else echo "Working braker_ab_initio folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/braker_ab_initio; date
fi

#--------------------------------------------------------------------------------

## Import master partition file 
guide_file=$WORKING_FOLDER_SCRATCH/scaffold_names_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 2, 18919; 631 partitions)
# Scaffold name       # Partition
# Backbone_10001              1
# Backbone_10003              1
# Backbone_10004              1
# Backbone_10005              1
# Backbone_10006              1
# ....

#--------------------------------------------------------------------------------


# Move to working directory
cd $WORKING_FOLDER_SCRATCH/braker_ab_initio

# Execute breaker in ab inition method (In this mode, GeneMark-ES is trained on the genome sequence, alone)
braker.pl \
--species=Nucella_canaliculata \
--genome=$REFERENCE \
--threads N 