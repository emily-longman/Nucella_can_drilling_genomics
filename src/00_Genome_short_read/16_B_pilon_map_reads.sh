#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Map_reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will map the reads with bwa mem.

# Load modules  
spack load gcc@9.3.0
spack load samtools@1.10

bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=$WORKING_FOLDER_SCRATCH/ntlink/final/final_assembly.ntLink.scaffolds.gap_fill.fa

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "pilon" ]
then echo "Working pilon folder exist"; echo "Let's move on."; date
else echo "Working pilon folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon; date
fi

# Move to pilon directory
cd $WORKING_FOLDER_SCRATCH/pilon

if [ -d "sams" ]
then echo "Working sams folder exist"; echo "Let's move on."; date
else echo "Working sams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/sams; date
fi

if [ -d "mapping_stats" ]
then echo "Working mapping_stats folder exist"; echo "Let's move on."; date
else echo "Working mapping_stats folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/mapping_stats; date
fi

if [ -d "bams" ]
then echo "Working bams folder exist"; echo "Let's move on."; date
else echo "Working bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/bams; date
fi

#--------------------------------------------------------------------------------

# Map reads to a reference

# This part will map reads to the reference genome. After reads have been mapped, they will be compressed into bam files, 
# sorted, and duplicates will be removed. I will also conduct an intermediary QC step with Qualimap. 
# Remember to take a look at the qualimap and the flagstat outputs to check for inconsistencies.

# Start pipeline

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/pilon

# Starting mapping
echo "Begin mapping" 
  
# I will conduct the mapping with BWA-MEM 2	
$bwa mem -M -t $CPU $REFERENCE \
$WORKING_FOLDER_SCRATCH/fastp/NC3_R1_clean.fastq.gz \
$WORKING_FOLDER_SCRATCH/fastp/NC3_R2_clean.fastq.gz \
> $WORKING_FOLDER_SCRATCH/pilon/sams/Ncan.sam

#--------------------------------------------------------------------------------

# I will now extract some summary stats
samtools flagstat --threads $CPU \
$WORKING_FOLDER_SCRATCH/pilon/sams/Ncan.sam \
> $WORKING_FOLDER_SCRATCH/pilon/mapping_stats/Ncan.flagstats_raw.sam.txt

# Build bam files
samtools view -b -q $QUAL --threads $CPU  \
$WORKING_FOLDER_SCRATCH/pilon/sams/Ncan.sam \
> $WORKING_FOLDER_SCRATCH/pilon/bams/Ncan.bam

#--------------------------------------------------------------------------------

# Inform that sample is done

echo "pipeline completed" $(date)