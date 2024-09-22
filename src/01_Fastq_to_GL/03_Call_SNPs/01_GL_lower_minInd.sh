#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Genotype_likelihoods_lower_minInd

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=02-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=200G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate genotype likelihoods of SNPs for all individuals

#--------------------------------------------------------------------------------

#Load modules 
spack load angsd@0.933
spack load samtools@1.10

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Path to the directory with the lane merged bams (filtered, sorted and duplicates removed). 
BAMS_FOLDER=$WORKING_FOLDER/bams_merged

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=40 #change accordingly in SLURM header
echo "using #CPUs ==" $NB_CPU

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods_all" ]
then echo "Working genotype_likelihoods_all folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods_all folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods_all; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/genotype_likelihoods_all

#--------------------------------------------------------------------------------

# Prepare bamlist
# This is a file with the name and full path of all the bam files to be processed.

# Move to bams folder
cd $BAMS_FOLDER

# Create bamlist for all Nucella samples
ls -d "$PWD/"* > $OUTPUT/Nucella_bam_lower_minInd.list

#--------------------------------------------------------------------------------

# Start pipeline

# Estimate Genotype Likelihoods's and allele frequencies for only the polymorphic sites. 
# Specifically calculate the SAF, MAF and GL for all individuals listed in the bam list.

# Move back to working directory
cd $WORKING_FOLDER

# File suffix to distinguish analysis choices
SUFFIX_2="SNPs_all_lower_minInd"

## Filter changes:
# minMapQ: threshold for minimum read mapping quality (Phred): increase from 20 to 30
# minInd: set min number of individuals to keep a site to 85% of 192 individuals
# setMinDepthInd: discard individual if sequencing depth for an individual is below 0.1
# skipTriallelic: don’t use sites with >2 alleles
# minMaf: Keep only sites with minor allele freq > some proportion (0.01)

# Generate GL's for polymorphic sites for all Nucella samples
angsd -b ${OUTPUT}/Nucella_bam_lower_minInd.list \
-ref ${REFERENCE} -anc ${REFERENCE} \
-out ${OUTPUT}/Nucella_${SUFFIX_2} \
-P $NB_CPU \
-doMaf 1 -doSaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-remove_bads 1 -baq 1 -skipTriallelic 1 -minMapQ 30 -minQ 20 \
-minInd 96 -setMinDepthInd 0.1 -minMaf 0.01 -setMaxDepth 600 \
-SNP_pval 1e-6 

# Change setMaxDepth to 600 (i.e., 3 * expected coverage (1X) * ~200 ind  )
# Switched GL to 2 rather than 1 (i.e., GL from GATK model rather than Samtools)

# -P: number of threads
# -doMaf 1: estimate allele frequencies
# -doSaf 1: estimate the SFS and/or neutrality tests genotype calling
# -GL 1: estimate genotype likelihoods (GL) using the Samtools formula (1)
# -doGlf 2: gives us the Beagle format which will be used by pcangsd
# -doMajorMinor 1: infer the major/minor using different approaches
# -doCounts 1: calculate various counts statistics

# -remove_bads 1: remove reads flagged as ‘bad’ by samtools
# -baq 1: estimates base alignment qualities for bases around indels
# -skipTriallelic 1: don’t use sites with >2 alleles
# -minMapQ 30: threshold for minimum read mapping quality (Phred)
# -minQ 20: threshold for minimum base quality (Phred)
# -C 50: enforce downgrading of map quality if contains excessive mismatches

# -minInd 163: min number of individuals to keep a site (~85%)
# -setMinDepthInd 0.1: min read depth for an individual to count towards a site
# -minMaf 0.01: Keep only sites with minor allele freq > some proportion.
# -setMaxDepth: Keep SNPs with a maximum total depth (typically set at 2-4 times the expected coverage times the number of ind to remove repeat regions)

# -SNP_pval 1e-6: Keep only site highly likely to be polymorphic (SNPs)
