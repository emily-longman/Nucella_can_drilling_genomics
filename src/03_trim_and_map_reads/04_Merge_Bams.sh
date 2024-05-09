#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Merge_bams

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
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will gather all sample data across the many lanes of sequencing

#Load modules 
module load samtools-1.10-gcc-7.3.0-pdbkohx
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#Folder for joint bams
JOINT_BAMS=/gpfs2/scratch/elongman/scratch/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/joint_bams

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=Merge_Bams

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

SAMPLE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/GuideFile.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##               File1                             File2              Snail_ID  Sample#  Lane#    Merged_name    Merged_bam_name
## FB1-1_S84_L002_R1_001.fastq.gz    FB1-1_S84_L002_R2_001.fastq.gz    FB1-1     S84     L002    FB1-1_S84_L002    FB1-1_S84
## FB1-1_S84_L007_R1_001.fastq.gz    FB1-1_S84_L007_R2_001.fastq.gz    FB1-1     S84     L007    FB1-1_S84_L007    FB1-1_S84
## FB1-1_S84_L008_R1_001.fastq.gz    FB1-1_S84_L008_R2_001.fastq.gz    FB1-1     S84     L008    FB1-1_S84_L008    FB1-1_S84
## FB1-2_S173_L002_R1_001.fastq.gz   FB1-2_S173_L002_R2_001.fastq.gz   FB1-2     S173    L002    FB1-2_S173_L002   FB1-2_S173
## ...
## MP9-10_S26_L007_R1_001.fastq.gz   MP9-10_S26_L007_R2_001.fastq.gz   MP9-10    S26     L007    MP9-10_S26_L007   MP9-10_S26
## MP9-10_S26_L008_R1_001.fastq.gz   MP9-10_S26_L008_R2_001.fastq.gz   MP9-10    S26     L008    MP9-10_S26_L008   MP9-10_S26

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

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "joint_bams" ]
then echo "Working joint_bams folder exist"; echo "Let's move on."; date
else echo "Working joint_bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/joint_bams; date
fi

if [ -d "joint_bams_qualimap" ]
then echo "Working joint_bams_qualimap folder exist"; echo "Let's move on."; date
else echo "Working joint_bams_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/joint_bams_qualimap; date
fi

if [ -d "Merged_Bams" ]
then echo "Working Merged_Bams folder exist"; echo "Let's move on."; date
else echo "Working Merged_Bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Merged_Bams; date
fi

if [ -d "Merged_Bams_qualimap" ]
then echo "Working Merged_Bams_qualimap folder exist"; echo "Let's move on."; date
else echo "Working Merged_Bams_qualimap folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Merged_Bams_qualimap; date
fi

#--------------------------------------------------------------------------------
# Start pipeline

# Here I will merge the bam outputs from the merge and unmerged portions of the pipeline. These will be named 'joint bams'
# Subsequently, I will once again sort and remove duplicated, before performing the final QC on the aligment.

# Merge bams
java -Xmx$JAVAMEM -jar $PICARD MergeSamFiles \
I=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam \
I=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam \
O=$WORKING_FOLDER/joint_bams/${i}.joint.bam

# Sort merge bams
java -Xmx$JAVAMEM -jar $PICARD SortSam \
I=$WORKING_FOLDER/joint_bams/${i}.joint.bam \
O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

# Remove duplicates of final file
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
I=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \
O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam  \
M=$WORKING_FOLDER/mapping_stats/${i}.joint.dupstat.txt \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=true

# Assess quality of final file
$qualimap bamqc \
-bam $WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam  \
-outdir $WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${i} \
--java-mem-size=$JAVAMEM
 
# Remove intermediary files
rm $WORKING_FOLDER/joint_bams/${i}.joint.bam
rm $WORKING_FOLDER/joint_bams/${i}.joint.srt.bam

#--------------------------------------------------------------------------------

# Here I will merge the joint bam outputs for the multiple lanes of sequencing. These will be named 'Lanes merged'

echo "I will merge these files" $JOINT_BAMS/${i}_*.joint.srt.rmdp.bam

#Make temporary linefile with list of input BAM files
ls $JOINT_BAMS/${i}_*.joint.srt.rmdp.bam > ${i}.guide.txt

samtools merge \
-b ${i}.guide.txt \
$WORKING_FOLDER/Merged_Bams/${i}.Lanes_merged.bam

#remove the temporary guide file
rm ${i}.guide.txt

# Assess quality of final file
$qualimap bamqc \
-bam $WORKING_FOLDER/Merged_Bams/${i}.Lanes_merged.bam \
-outdir $WORKING_FOLDER/Merged_Bams_qualimap/Qualimap_LaneMerged_${i} \
--java-mem-size=$JAVAMEM

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline completed" $(date)