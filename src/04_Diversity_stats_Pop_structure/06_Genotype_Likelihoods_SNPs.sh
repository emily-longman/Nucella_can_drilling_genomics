#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Genotype_likelihoods_SNPs

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=06-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=500G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------


#Load modules 
module load angsd-0.933-gcc-7.3.0-4wsdzjw
module load samtools-1.10-gcc-7.3.0-pdbkohx

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#Folder for joint bams
BAMS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/Merged_Bams

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#Name of pipeline
PIPELINE=Genotype_Likelihoods_SNPs

#--------------------------------------------------------------------------------
# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#--------------------------------------------------------------------------------
# Move to working directory
cd $WORKING_FOLDER

# Begin Pipeline

# This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "Your unique run id is:" $unique_run_id

echo $PIPELINE
echo $WORKING_FOLDER

if [[ -e "${PIPELINE}.warnings.log" ]]
then echo "Warning log exist"; echo "Let's move on."; date
else echo "Warning log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/${PIPELINE}.warnings.log; date
fi

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genotype_likelihoods_SNPs" ]
then echo "Working genotype_likelihoods_SNPs folder exist"; echo "Let's move on."; date
else echo "Working genotype_likelihoods_SNPs folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/genotype_likelihoods_SNPs; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/genotype_likelihoods_SNPs

#--------------------------------------------------------------------------------

## PREPARE bamlist
# This is a file with the name and full path of all the bam files to be processed.

# Move to bams folder
cd $BAMS_FOLDER

# Create bamlist for all Nucella samples
ls -d "$PWD/"* > $OUTPUT/Nucella_bam.list

#--------------------------------------------------------------------------------

# Estimate Genotype Likelihoods's and allele frequencies for only the polymorphic sites

# Move back to working directory
cd $WORKING_FOLDER

# File suffix to distinguish analysis choices
SUFFIX_2="SNPs"

# Generate GL's for polymorphic sites for all Nucella samples
angsd -b ${OUTPUT}/Nucella_bam.list \
-ref ${REFERENCE} \
-anc ${REFERENCE} \
-out ${OUTPUT}/Nucella_${SUFFIX_2} \
-nThreads $CPU \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-GL 1 \
-doSaf 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doGlf 2 \
-SNP_pval 1e-6 \
-minMaf 0.01


# -doSaf 1	output allele frequency likelihoods for each site
# -doCounts 1	output allele counts for each site
# -doMajorMinor 1	fix major and minor alleles the same across all samples
# -doMaf 1	calculate minor allele frequencies
# -SNP_pval 1e-6	Keep only site highly likely to be polymorphic (SNPs)
# -minMaf 0.01	Keep only sites with minor allele freq > some proportion.

# -doGlf 2 gives us the Beagle format which will be used by pcangsd