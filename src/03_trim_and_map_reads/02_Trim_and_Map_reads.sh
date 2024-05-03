#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Trim_and_Map

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

# This script is the second step in the trim and map pipeline

#Load modules 
module load gcc/10.5.0
module load bwa-0.7.17-gcc-7.3.0-terdbma
module load fastqc-0.11.7-gcc-7.3.0-vcaesw7
module load samtools-1.10-gcc-7.3.0-pdbkohx

bbduk=/gpfs1/home/e/l/elongman/software/bbmap/bbduk.sh 
qualimap=/gpfs1/home/e/l/elongman/software/qualimap_v2.3/qualimap
PICARD=/gpfs1/home/e/l/elongman/software/picard.jar

#Define important file locations
#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#This is a unique number id which identifies this run
unique_run_id=`date +%N`

#Name of pipeline
PIPELINE=Trim_and_Map

#--------------------------------------------------------------------------------
#Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18g # Java memory

#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

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
i=`awk -F "\t" '{print $6}' $SAMPLE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo $i

#--------------------------------------------------------------------------------

# Index reference (this step only needs to be done once)

# This indexing step only needs to be done once for the reference file.
#cd /netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024
#bwa index $REFERENCE 
# -p is the preix of the output databse 
# -a is the algorithm for constructing BWT index ('is' - linear-time algorithm for constructing suffix array; bwtsw 

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Begin Pipeline

#This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "Your unique run id is:" $unique_run_id

echo $PIPELINE
echo $WORKING_FOLDER

if [[ -e "${PIPELINE}.warnings.log" ]]
then echo "Warning log exist"; echo "Let's move on"; date
else echo "Warning log doesnt exist. Let's fix that"; touch $WORKING_FOLDER/${PIPELINE}.warnings.log; date
fi

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on"; date
else echo "Completion log doesnt exist. Let's fix that"; touch $WORKING_FOLDER/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "mapping_stats" ]
then echo "Working mapping_stats folder exist"; echo "Let's move on"; date
else echo "Working mapping_stats folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/mapping_stats; date
fi

if [ -d "read_stats" ]
then echo "Working read_stats folder exist"; echo "Let's move on"; date
else echo "Working read_stats folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/read_stats; date
fi

if [ -d "joint_bams" ]
then echo "Working joint_bams folder exist"; echo "Let's move on"; date
else echo "Working joint_bams folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/joint_bams; date
fi

if [ -d "joint_bams_qualimap" ]
then echo "Working joint_bams_qualimap folder exist"; echo "Let's move on"; date
else echo "Working joint_bams_qualimap folder doesnt exist. Let's fix that"; mkdir $WORKING_FOLDER/joint_bams_qualimap; date
fi

#--------------------------------------------------------------------------------
# Start pipeline
# Lets do some light trimming of the reads


### Work on merged reads
### Decide on the trimming parameters based on fastQC step done before this script.

echo ${i} "Trimming merged reads"

$bbduk \
in=`echo $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq` \
out=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.trim.fq \
ftl=12 ftr=288 qtrim=w trimq=20

rm  $WORKING_FOLDER/merged_reads/${i}/${i}.merged.reads.strict.fq

### Work on unmerged reads
echo ${i} "Trimming unmerged reads"
	
$bbduk \
in=`echo $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq` \
in2=`echo $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq` \
out=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.1.fq \
out2=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.2.fq \
ftl=12 qtrim=w trimq=20
	
rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.1.fq
rm $WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.2.fq


#--------------------------------------------------------------------------------
# Map reads to a reference

# This part will map reads to the reference genome. Because the reads are likely split into two groups, 
# this script will loop over both types of reads. After reads have been mapped, they will be compressed into bam files, 
# sorted, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap. 
# Because there are inherent QC steps here, I have avoided adding extra "warnings" in the log. 
# Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.  

for j in merged unmerged
do # Begin loop of j

########################################
#J loop#	# Starting mapping
echo "I will first map ${j} reads of" ${i}
	
#J loop#	# I will conduct the mapping with BWA-MEM	

if [[ ${j} == "merged" ]]; 
then echo "seems this is merged data, lets map it"; 
bwa mem -M -t $CPU $REFERENCE \ 
$WORKING_FOLDER/merged_reads/${i}/${i}.${j}.reads.strict.trim.fq > $WORKING_FOLDER/merged_reads/${i}/${i}.${j}.sam

elif [[ ${j} == "unmerged" ]]; 
then echo "seems this is unmerged data, lets map it using a 1-2 approach"; 
bwa mem -M -t $CPU $REFERENCE \ 
$WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.reads.trim.1.fq \
$WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.reads.trim.2.fq \
> $WORKING_FOLDER/unmerged_reads/${i}/${i}.${j}.sam

$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.reads.trim.1.fq

else echo "I cant tell what type of data this is -- WARNING!"; echo ${i} "Something is wrong at the mapping stage" $(date) \ 
$Project_name.warnings.$unique_run_id.log
fi

#J loop#	#I will now extract some summary stats
samtools flagstat --threads $CPU \ 
$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
> $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt

#J loop#	#build bam files
samtools view -b -q $QUAL --threads $CPU  \ 
$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
> $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam

done # End loop of j

# Inform that sample is done
###########################################################################
###########################################################################
# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline completed" $(date)