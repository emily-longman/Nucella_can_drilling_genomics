#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=pilon_polish

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=6-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=800G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Call package (installed with conda)
spack load samtools@1.10
PILONJAR=/gpfs1/home/e/l/elongman/software/pilon-1.24.jar

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=$WORKING_FOLDER_SCRATCH/ntlink/final/final_assembly.ntLink.scaffolds.gap_fill.fa

# This is the location of the cleaned and indexed bams
BAM=$WORKING_FOLDER_SCRATCH/pilon/bams_clean/Ncan.srt.rmdp.bam

#--------------------------------------------------------------------------------

# Java memory limits
JAVA_TOOL_OPTIONS="-Xmx300G -Xss2560k" # set maximum heap size to something reasonable

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/pilon

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "genome_chunks" ]
then echo "Working genome_chunks folder exist"; echo "Let's move on."; date
else echo "Working genome_chunks folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/pilon/genome_chunks; date
fi

#--------------------------------------------------------------------------------

# Use pilon to polish the genome 

java -Xmx800G -jar $PILONJAR --genome $REFERENCE --frags $BAM --output N.canaliculata_polished_genome --outdir $WORKING_FOLDER_SCRATCH/pilon/polished_genome

# --frags for paired-end sequencing of DNA fragments, such as Illumina paired-end reads of fragment size <1000bp.

#--------------------------------------------------------------------------------

conda deactivate