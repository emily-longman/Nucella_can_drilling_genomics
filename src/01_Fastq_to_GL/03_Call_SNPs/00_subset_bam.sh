#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Subset_bams_GL_test

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
#SBATCH --array=1-192%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

#Load modules 
spack load samtools@1.10
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Path to the directory with the lane merged bams (filtered, sorted and duplicates removed). 
#Note I copied all of the bams from the original bam directory ("bams_merged") into "bams_merged_indexed"
BAMS_FOLDER=$WORKING_FOLDER/bams_merged_indexed

#--------------------------------------------------------------------------------

# Change to working folder
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "bams_merged_subset" ]
then echo "Working bams_merged_subset folder exist"; echo "Let's move on."; date
else echo "Working bams_merged_subset folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/bams_merged_subset; date
fi

#--------------------------------------------------------------------------------
## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

GUIDE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/Guide_Files/Guide_File_merge.txt

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
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo $i

#--------------------------------------------------------------------------------

# Index Bam files

java -jar $PICARD BuildBamIndex \
I=$WORKING_FOLDER/bams_merged_indexed/${i}.Lanes_merged.bam \
O=$WORKING_FOLDER/bams_merged_indexed/${i}.Lanes_merged.bam.bai

#samtools index $WORKING_FOLDER/bams_merged_indexed/${i}.Lanes_merged.bam

#--------------------------------------------------------------------------------

# Samtools coverage to determine scaffold coverage

#Make temporary linefile with list of input BAM files
#ls $WORKING_FOLDER/bams_merged/*.Lanes_merged.bam > $WORKING_FOLDER/Guide_Files/bams_cov_guide.txt

#samtools coverage -m -b $WORKING_FOLDER/Guide_Files/bams_cov_guide.txt 

#--------------------------------------------------------------------------------

# Testing 

# Subset the data for only a few scaffolds
# This code will make a bam file for each backbone

SCAFFOLD_LIST="Backbone_537 Backbone_552 Backbone_656 Backbone_800 Backbone_980 Backbone_985"

for s in $SCAFFOLD_LIST
do 
echo "Processing $s"
samtools view -b $WORKING_FOLDER/bams_merged_indexed/${i}.Lanes_merged.bam ${s} > $WORKING_FOLDER/bams_merged_subset/${i}.${s}.subset.bam
done

# If only want one scaffold use the code below:
#samtools view -b $WORKING_FOLDER/bams_merged_indexed/${i}.RG.bam "ntLink_33941" > $WORKING_FOLDER/bams_merged_indexed/${i}.subset.bam

#--------------------------------------------------------------------------------

# Merge the bam files for each backbone together, so that there is only one bam file for each sample 

#Make temporary linefile with list of input BAM files
ls $WORKING_FOLDER/bams_merged_subset/${i}.*.subset.bam > ${i}.guide.txt

samtools merge \
-b ${i}.guide.txt \
$WORKING_FOLDER/bams_merged_subset/${i}.subset.bam

#Remove the temporary guide file
rm ${i}.guide.txt

#Remove the individual scaffold files
cd $WORKING_FOLDER/bams_merged_subset
rm ${i}.Backbone_*.subset.bam