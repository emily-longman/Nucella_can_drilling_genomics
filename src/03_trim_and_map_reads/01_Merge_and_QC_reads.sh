#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Merge_and_QC_reads

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
#SBATCH --array=1

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will initiate a pipeline which will do some quality QC on the reads and then will proceed to map the reads to a reference genome.

# Notes on Nomenclature:
# The first part of the script split raw reads into insert-"merged"-reads (hereby called merged) and unmerged reads (those which were not merged). 
# As such, all operations done using either of these reads will have the term "merged" or "unmerged" attached to them. 
# At a later point in the script, I combine bam files using "samtools merge" the output of this combination is a joint-bam file (hereby called "joint"). 
# Thus, the joint suffix referes to this step. Other suffix used here are: "srt" which mean picard sorted, and "rmdp" which mean picard-removed duplicated reads.
  
#Load modules 
module load gcc/10.5.0
module load bwa-0.7.17-gcc-7.3.0-terdbma
module load fastqc-0.11.7-gcc-7.3.0-vcaesw7
module load samtools-1.10-gcc-7.3.0-pdbkohx

bbmerge=/gpfs1/home/e/l/elongman/software/bbmap/bbmerge.sh #executable

#Define important file locations
#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/FL_2000/Assembly.fasta

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=Merge_and_Trim

#--------------------------------------------------------------------------------
#Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18g # Java memory

#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. one sample name per line with all the infor
SAMPLE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/GuideFile.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
##               File1                             File2              Snail_ID  Sample#  Lane#    Merged_name
## FB1-1_S84_L002_R1_001.fastq.gz    FB1-1_S84_L002_R2_001.fastq.gz    FB1-1     S84    L002    FB1-1_S84_L002 
## FB1-1_S84_L007_R1_001.fastq.gz    FB1-1_S84_L007_R2_001.fastq.gz    FB1-1     S84    L007    FB1-1_S84_L007
## FB1-1_S84_L008_R1_001.fastq.gz    FB1-1_S84_L008_R2_001.fastq.gz    FB1-1     S84    L008    FB1-1_S84_L008
## FB1-2_S173_L002_R1_001.fastq.gz   FB1-2_S173_L002_R2_001.fastq.gz   FB1-2     S173   L002    FB1-2_S173_L002
## ...
## MP9-10_S26_L007_R1_001.fastq.gz   MP9-10_S26_L007_R2_001.fastq.gz   MP9-10    S26    L007    MP9-10_S26_L007
## MP9-10_S26_L008_R1_001.fastq.gz   MP9-10_S26_L008_R2_001.fastq.gz   MP9-10    S26    L008    MP9-10_S26_L008


#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
 
i=`awk -F "\t" '{print $6}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`
read1=`awk -F "\t" '{print $1}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`
read2=`awk -F "\t" '{print $2}' $SAMPLE_FILE | sed -n ${SLURM_ARRAY_TASK_ID}p`

#--------------------------------------------------------------------------------
# Begin Pipeline

#This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "Your unique run id is:" $unique_run_id

if [[ -e "${PIPELINE}.warnings.log" ]]
then
	echo "Warning log exist"
	echo "Lets move on"
	date
else 
	echo "Warning log doesnt exist. Lets fix that"
	touch $WORKING_FOLDER/${PIPELINE}.warnings.log
	date
fi

if [[ -e "${PIPELINE}.completion.log" ]]
then
	echo "Completion log exist"
	echo "Lets move on"
	date
else 
	echo "Completion log doesnt exist. Lets fix that"
	touch $WORKING_FOLDER/${PIPELINE}.completion.log
	date
fi

#--------------------------------------------------------------------------------
# Generate Folders and Files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

# Generating new folders 
echo "have you checked if the folders were already built with mkdir?"
if [[ -d "merged_reads" ]]
then
	echo "Working merged_reads folder exist"
	echo "Lets move on"
	date
else 
	echo "Working merged_reads folder doesnt exist. Lets fix that"
	mkdir $WORKING_FOLDER/merged_reads
	date
fi

if [ -d "unmerged_reads" ]
then
	echo "Working unmerged_reads folder exist"
	echo "Lets move on"
	date
else 
	echo "Working unmerged_reads folder doesnt exist. Lets fix that"
	mkdir $WORKING_FOLDER/unmerged_reads
	date
fi

if [ -d "fastqc_merged" ]
then
	echo "Working fastqc_merged folder exist"
	echo "Lets move on"
	date
else 
	echo "Working fastqc_merged folder doesnt exist. Lets fix that"
	mkdir $WORKING_FOLDER/fastqc_merged
	date
fi

#--------------------------------------------------------------------------------
#  Merge and QC reads

# This part of the pipeline will merge the reads. It is very likely that the reads will be split into merged and unmerged. 
# Both reads will be mapped. This loop operates using a while-read-do-done structure. The while loop is feed a file "SAMPLE_FILE", 
# where  all sample names are stored, one name per line. This can be leveraged for parallelization.

# Move to working directory
cd $WORKING_FOLDER

	echo ${i} "is now processing"
	date

	mkdir $WORKING_FOLDER/merged_reads/${i}
	mkdir $WORKING_FOLDER/unmerged_reads/${i}
	
	echo "now merging reads for" ${i}
	
	$bbmerge \
	in1=${RAW_READS}/$read1 in2=${RAW_READS}/$read2 \
	out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq \
	outu1=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq \
	outu2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq \
	-strict
	
		#Sanity checks	
	if [ -s $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq ]; then
	echo ${i} "merged reads file is not empty... thats good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Merged reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi
	
	if [ -s $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq ]; then
	echo ${i} "Pair 1 reads file is not empty... thats good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Pair 1 reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi
	
	if [ -s $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq ]; then
	echo ${i} "Pair 2 reads file is not empty... thats good"
	else
	echo "File is empty -- WARNING ISSUED!"
	echo ${i} "Pair 2 reads is empty! check the script, debug, and rerun" >> $WORKING_FOLDER/${PIPELINE}.warnings.log
	fi	

#--------------------------------------------------------------------------------
# Now lest do some QC on the reads

	fastqc $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq \
	--outdir $WORKING_FOLDER/fastqc_merged


#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will produce a notification stating the completion of the script. 

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline" ${PIPELINE} $(date)