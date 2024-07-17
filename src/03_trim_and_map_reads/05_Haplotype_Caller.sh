#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Haplotype_Caller

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Submit job array
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/HaploCaller.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will initiate a pipeline which will add read group info and index bams. It will then proceed to call haplotypes (gVCFs)

#Load modules 
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
#gatk=/gpfs1/home/e/l/elongman/software/gatk-4.5.0.0/gatk # VACC's verion of java isn't up to date enough to be compatible
gatk=/gpfs1/home/e/l/elongman/software/gatk-4.2.6.0/gatk
spack load htslib@1.10.2 #tabix
spack load samtools@1.10

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#Folder for joint bams
BAMS_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/Merged_Bams

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=Haplocaller

#Sample suffixes and post-fixes. What tags are expected across all samples?
# This comes from the previous pipeline
SUFFIX="Lanes_merged"

#--------------------------------------------------------------------------------
# Define parameters

# Java parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

#Read Information
Group_library="Longman_2023"
Library_Platform="illumina"
Group_platform="L2023"

#HaploCaller -- heterozygocity prior
HET=0.005

#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

SAMPLE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/GuideFileMerged.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
## Snail_ID  Sample#   Merged_name 1    Merged_name 2    Merged_name 3   Merged_bam_name
##  FB1-1     S84     FB1-1_S84_L002   FB1-1_S84_L007   FB1-1_S84_L008     FB1-1_S84
##  FB1-2     S173    FB1-2_S173_L002  FB1-2_S173_L007  FB1-2_S173_L008    FB1-2_S173
##  FB1-5     S109    FB1-5_S109_L002  FB1-5_S109_L007  FB1-5_S109_L008    FB1-5_S109
##  ...
##  MP9-9     S191    MP9-9_S191_L002  MP9-9_S191_L007  MP9-9_S191_L008    MP9-9_S191
##  MP9-10    S26     MP9-10_S26_L002  MP9-10_S26_L007  MP9-10_S26_L008    MP9-10_S26

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $6}' $SAMPLE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo $i

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

#--------------------------------------------------------------------------------

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "RGSM_final_bams" ]
then echo "Working RGSM_final_bams folder exist"; echo "Let's move on."; date
else echo "Working RGSM_final_bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/RGSM_final_bams; date
fi

if [ -d "haplotype_calling" ]
then echo "Working haplotype_calling folder exist"; echo "Let's move on."; date
else echo "Working haplotype_calling folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/haplotype_calling; date
fi

#--------------------------------------------------------------------------------

# Forcing a uniform read group to the joint bam file

java -jar $PICARD AddOrReplaceReadGroups \
I=$BAMS_FOLDER/${i}.$SUFFIX.bam \
O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
RGLB=$Group_library \
RGPL=$Library_Platform \
RGPU=$Group_platform \
RGSM=${i}

#--------------------------------------------------------------------------------

# Index Bam files

java -jar $PICARD BuildBamIndex \
I=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam.bai

#samtools index $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam

#--------------------------------------------------------------------------------

# Haplotype Calling

# Need .dict and .fai index files generated from Picard and samtools (in '02_Index_reference' step)

# Call haplotypes with GATK

$gatk --java-options "-Xmx${JAVAMEM}" HaplotypeCaller \
-R $REFERENCE \
-I $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
-O $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf \
--heterozygosity $HET \
-ERC GVCF

#--------------------------------------------------------------------------------

# Compress and index with Tabix (within HTSlib in Samtools)

bgzip $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf
tabix $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf.gz

echo ${i} "completed" $(date) >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "done" $(date)