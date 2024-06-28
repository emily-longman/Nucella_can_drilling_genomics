#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Subset_bams

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Submit job array
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

#Load modules 
module load samtools-1.10-gcc-7.3.0-pdbkohx
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=Subset_bams

#--------------------------------------------------------------------------------
# Define parameters

# Java parameters
#CPU=$SLURM_CPUS_ON_NODE
#echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
#QUAL=40 # Quality threshold for samtools
#JAVAMEM=18G # Java memory

#Read Information
Group_library="Longman_2023"
Library_Platform="illumina"
Group_platform="L2023"

#--------------------------------------------------------------------------------
# Begin Pipeline

# This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "Your unique run id is:" $unique_run_id

echo $PIPELINE
echo $WORKING_FOLDER

cd $WORKING_FOLDER/Logs

# Make Warning logs
if [[ -e "${PIPELINE}.warnings.log" ]]
then echo "Warning log exist"; echo "Let's move on."; date
else echo "Warning log doesnt exist. Let's fix that."; touch $$WORKING_FOLDER/Logs/${PIPELINE}.warnings.log; date
fi

#--------------------------------------------------------------------------------

# Change to working folder
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "RGSM_final_bams" ]
then echo "Working RGSM_final_bams folder exist"; echo "Let's move on."; date
else echo "Working RGSM_final_bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/RGSM_final_bams; date
fi

if [ -d "BAMS_subset" ]
then echo "Working BAMS_subset folder exist"; echo "Let's move on."; date
else echo "Working BAMS_subset folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/BAMS_subset; date
fi

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

# Change to output folder
#cd $WORKING_FOLDER/BAMS_subset

# Samtools coverage to determine scaffold coverage

#Make temporary linefile with list of input BAM files
#ls $WORKING_FOLDER/Merged_Bams/*.Lanes_merged.bam > guide.txt

#samtools coverage -m -b guide.txt

#--------------------------------------------------------------------------------

# Examples form online
#samtools view -b input.bam "Chr10:18000-45500" > output.bam

# Index bam files

# Testing 
#samtools view -b $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam "ntLink_33941" > $WORKING_FOLDER/BAMS_subset/${i}.subset.bam

#Can't figure out how to specify multiple ntLink scaffolds