#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Get_sites

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=4:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will extract the sites and region files for MAF scripts
# 1) extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
# 2) order the sites file by scaffold names 
# 3) makes a region file matching the sites files and with same order
# 4) index sites file

# Load modules
spack load angsd@0.933

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Output folder
GL_FOLDER=$WORKING_FOLDER/genotype_likelihoods_all 

#--------------------------------------------------------------------------------


echo "from the maf file, extract a list of SNP chr, position, major all, minor all"
gunzip $GL_FOLDER/Nucella_SNPs_all.mafs.gz

INFILE=$GL_FOLDER/Nucella_SNPs_all.mafs

#THESE ARENT UPDATED:
OUTFILE_sites=02_info/sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"
OUTFILE_regions=02_info/regions_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR"

Rscript Rscripts/make_sites_list_maxdepth_simple.R "$INFILE" "$OUTFILE_sites" "$OUTFILE_regions"

angsd sites index $OUTFILE_sites
