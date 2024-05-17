#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Check_Coverage #Accidentally still had "Merge & QC"

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

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
module load gcc/10.5.0
module load bwa-0.7.17-gcc-7.3.0-terdbma
module load fastqc-0.11.7-gcc-7.3.0-vcaesw7
module load samtools-1.10-gcc-7.3.0-pdbkohx
bbduk=/gpfs1/home/e/l/elongman/software/bbmap/bbduk.sh 
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

#Define important file locations
#RAW READS indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#--------------------------------------------------------------------------------
#Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
QUAL=40 # Quality threshold for samtools
JAVAMEM=18g # Java memory

#--------------------------------------------------------------------------------

WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Check_Coverage

# Cat all reads

# Concatenate reads into one fastq.gz file
cat $RAW_READS/*.fastq.gz > $WORKING_FOLDER/all.illumina.fastq.gz

# Map reads
$bwa mem -t $CPU $REFERENCE \
$WORKING_FOLDER/all.illumina.fastq.gz \
> $WORKING_FOLDER/all.illumina.sam

# Summary stats
samtools flagstat --threads $CPU \
$WORKING_FOLDER/all.illumina.sam \
> $WORKING_FOLDER/flagstats.all.illumina.txt

# Build bam files
samtools view -b --threads $CPU  \
$WORKING_FOLDER/all.illumina.sam \
> $WORKING_FOLDER/all.illumina.bam


# Sort with picard
# Notice that once a file has been sorted it is added the "srt" suffix
java -Xmx$JAVAMEM -jar $PICARD SortSam \
I=$WORKING_FOLDER/all.illumina.bam \
O=$WORKING_FOLDER/all.illumina.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

# Lets do QC on the bam file
$qualimap bamqc \
-bam $WORKING_FOLDER/all.illumina.srt.bam \
-outdir $WORKING_FOLDER/Qualimap.all.illumina.srt \
--java-mem-size=$JAVAMEM