#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Filter_bams

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Submit job array
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Filter_Bams.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

#Load modules 
spack load samtools@1.10
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=Filter_Bams

#--------------------------------------------------------------------------------
# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

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

cd Logs

if [[ -e "${PIPELINE}.warnings.log" ]]
then echo "Warning log exist"; echo "Let's move on."; date
else echo "Warning log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/Logs/${PIPELINE}.warnings.log; date
fi

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/Logs/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "Filtered_Bams" ]
then echo "Working Filtered_Bams folder exist"; echo "Let's move on."; date
else echo "Working Filtered_Bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Filtered_Bams; date
fi

if [ -d "Filtered_Bams_qualimap" ]
then echo "Working Filtered_Bams_qualimap folder exist"; echo "Let's move on."; date
else echo "Working Filtered_Bams_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Filtered_Bams_qualimap; date
fi

#--------------------------------------------------------------------------------

# Filter merged bam files with samtools view and add flags
# -q = Skip alignments with MAPQ smaller than $QUAL (40)
# 0x0002 = read mapped in proper pair (0x2)*
# 0x0004 = read unmapped (0x4)
# 0x0008 = mate unmapped (0x8)*

samtools view \
-b \
-q $QUAL \
-f 0x0002 -F 0x0004 -F 0x0008 \
--threads $CPU  \
$WORKING_FOLDER/Merged_Bams/${i}.Lanes_merged.bam > $WORKING_FOLDER/Filtered_Bams/${i}.flt.bam


# Assess quality of final file
$qualimap bamqc \
-bam $WORKING_FOLDER/Filtered_Bams/${i}.flt.bam \
-outdir $WORKING_FOLDER/Merged_Bams_qualimap/Qualimap_LaneMerged_${i} \
--java-mem-size=$JAVAMEM

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline completed" $(date)