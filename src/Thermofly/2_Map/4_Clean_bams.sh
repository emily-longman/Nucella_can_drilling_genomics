#!/usr/bin/env bash  
#  
#SBATCH -J Clean_bams  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 40G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-22
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Clean bams

# Load software  
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

#--------------------------------------------------------------------------------

echo ${SLURM_ARRAY_TASK_ID}

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.vNov11.2024.tsv
ref=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_nanopore.fasta.masked

#--------------------------------------------------------------------------------

# Use metadata file to extract sample names
SAMP_NAME=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo $SAMP_NAME

#--------------------------------------------------------------------------------

# Parameters for software
CPU=6
QUAL=40 
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir bams_clean

#--------------------------------------------------------------------------------

# Clean the bam files

# Filter bam files and add flags
samtools view \
-b \
-q $QUAL \
-f 0x0002 -F 0x0004 -F 0x0008 \
--threads $CPU  \
$working_folder/bams/${SAMP_NAME}.bam \
> $working_folder/bams_clean/${SAMP_NAME}.bam

# Sort with picard
java -Xmx$JAVAMEM -jar $PICARD SortSam \
I=$working_folder/bams_clean/${SAMP_NAME}.bam \
O=$working_folder/bams_clean/${SAMP_NAME}.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

# Remove duplicates with picard
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \
I=$working_folder/bams_clean/${SAMP_NAME}.srt.bam \
O=$working_folder/bams_clean/${SAMP_NAME}.srt.rmdp.bam \
M=$working_folder/bams_clean/${SAMP_NAME}.dupstat.txt \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

# Index with samtools
samtools index $working_folder/bams_clean/${SAMP_NAME}.srt.rmdp.bam

#--------------------------------------------------------------------------------

# Housekeeping and remove intermediate files
mv $working_folder/bams_clean/${SAMP_NAME}.dupstat.txt \
$working_folder/mapping_stats
rm $working_folder/bams_clean/${SAMP_NAME}.bam
rm $working_folder/bams_clean/${SAMP_NAME}.srt.bam

#--------------------------------------------------------------------------------

# Inform completion of pipeline
echo "pipeline completed" $(date)