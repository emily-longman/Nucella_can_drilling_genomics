#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=angsd_config

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss
#SBATCH --time=1:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_mask/N.canaliculata_assembly.fasta.masked

#Path to the directory with the lane merged bams (filtered, sorted and duplicates removed). 
BAMS_FOLDER=$WORKING_FOLDER/bams_merged

#--------------------------------------------------------------------------------

# Define parameters that will be used in the following scripts

# Filter : will keep SNP with minor allele freq above some proportion (over all individuals)
MIN_MAF=0.05

# Filter : will keep positions with at least MIN_DEPTH reads for each individual 
# This is not necessaily for all individuals, we consider a PERCENT_IND (percentage of individuals over all individuals in step 03, and within each pop at step 07)
# Advice: For min depth use a value that is a bit below what you expected. 
# Advice: For percent individual, avoid going below 50% and also consider the whole number of individuals. (it may makes sense to use 50% with 100 ind/pop, but you may want 90% with 9 ind/pop
PERCENT_IND=0.5
MIN_DEPTH=0.1

# Filter: will keep SNP with at least a coverage of this factor multiplied by the number of ind - across all ind. 
# Advice: we usually set it at 2-4 times the expected coverage to remove repeated regions
MAX_DEPTH_FACTOR=2

# Window size for sliding window FST & Thetas
WINDOW=25000 

# Window step
WINDOW_STEP=5000 

# Min number of pop to consider for NGS admix
K_MIN=1

# Maximum number of pop to consider for NGS admix
K_MAX=3