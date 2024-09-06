#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=index_genome

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

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

# This script with index the reference genome.

# Note: this step only needs to be done once

#Load modules 
spack load samtools@1.10
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

# Genome from scaffolding
ASSEMBLY=$WORKING_FOLDER_SCRATCH/ntlink/ntlink_k.30_w.200_r.6/final_assembly.fasta.k30.w200.z1000.ntLink.ntLink.ntLink.ntLink.ntLink.gap_fill.fa.k30.w200.z1000.ntLink.scaffolds.gap_fill.fa

#--------------------------------------------------------------------------------

# Move to the directory that the genome is currently stored
cd $WORKING_FOLDER_SCRATCH/ntlink/ntlink_k.30_w.200_r.6

# Index database sequences in the FASTA format 
$bwa index $ASSEMBLY 

# Generate the FASTA sequence dictionary file
java -jar $PICARD CreateSequenceDictionary \
R=$ASSEMBLY

# Generate the fasta index file for the reference
samtools faidx $ASSEMBLY