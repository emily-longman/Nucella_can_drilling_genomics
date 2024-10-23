#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Genotype_likelihoods

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

# This script will calculate genotype likelihoods of SNPs for all individuals.

#--------------------------------------------------------------------------------

#Load modules 
spack load angsd@0.933
spack load samtools@1.10

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_mask/N.canaliculata_assembly.fasta.masked

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

# Path to bam list.
BAM_LIST=$WORKING_FOLDER/guide_files/Nucella_bam.list

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=40 #change accordingly in SLURM header
echo "using #CPUs ==" $NB_CPU

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

# Extract parameters from config file
N_IND=$(wc -l $BAM_LIST | cut -d " " -f 1) 
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods_all" ]
then echo "Working genotype_likelihoods_all folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods_all folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods_all; date
fi

#--------------------------------------------------------------------------------

# Start pipeline

# Estimate Genotype Likelihoods's and allele frequencies for only the polymorphic sites. 
# Specifically calculate the SAF, MAF and GL for all individuals listed in the bam list.

# Move back to working directory
cd $WORKING_FOLDER

# Generate GL's for polymorphic sites for all Nucella samples
angsd \
-b $BAM_LIST \
-ref ${REFERENCE} -anc ${REFERENCE} \
-out $WORKING_FOLDER/genotype_likelihoods_all/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR" \
-P $NB_CPU \
-doMaf 1 -doSaf 1 -GL 2 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-remove_bads 1 -baq 1 -skipTriallelic 1 -uniqueOnly 1 -only_proper_pairs 1 -minMapQ 30 -minQ 20 -C 50 \
-minInd $MIN_IND -setMinDepthInd $MIN_DEPTH -minMaf $MIN_MAF -setMaxDepth $MAX_DEPTH \
-SNP_pval 1e-4 

# -P: number of threads

# -doMaf 1: estimate allele frequencies
# -doSaf 1: estimate the SFS and/or neutrality tests genotype calling
# -GL 2: estimate genotype likelihoods (GL) using the GATK formula
# -doGlf 2: gives us the Beagle format which will be used by pcangsd
# -doMajorMinor 1: infer the major/minor using different approaches
# -doCounts 1: calculate various counts statistics

# -remove_bads 1: remove reads flagged as ‘bad’ by samtools
# -baq 1: estimates base alignment qualities for bases around indels
# -skipTriallelic 1: don’t use sites with >2 alleles
# -uniqueOnly 1: Remove reads that have multiple best hits
# -only_proper_pairs 1: Include only proper pairs (pairs of read with both mates mapped correctly)
# -minMapQ 30: threshold for minimum read mapping quality (Phred)
# -minQ 20: threshold for minimum base quality (Phred)
# -C 50: enforce downgrading of map quality if contains excessive mismatches

# -minInd: min number of individuals to keep a site
# -setMinDepthInd: min read depth for an individual to count towards a site
# -minMaf: Keep only sites with minor allele freq > some proportion
# -setMaxDepth: Keep SNPs with a maximum total depth 

# -SNP_pval 1e-4: Keep only site highly likely to be polymorphic (SNPs)
