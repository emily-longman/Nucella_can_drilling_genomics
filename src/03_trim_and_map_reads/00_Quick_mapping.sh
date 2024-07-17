#!/usr/bin/env bash
#
#SBATCH -J dsubwa # A single job name for the array
#SBATCH -c 5
#SBATCH -N 1 # on one node
#SBATCH -t 20:00:00 #<= this may depend on your resources
#SBATCH --mem 60G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/hic.%A_%a.out # Standard output
#SBATCH -p bluemoon
#SBATCH --array=1-20

JAVAMEM=59G
CPU=5
SAMP=${SLURM_ARRAY_TASK_ID}

### programs
spack load samtools@1.10

### Link softwares
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
qualimap=/netfiles/nunezlab/Shared_Resources/Software/qualimap_v2.2.1/qualimap

###Link files
rloc=/netfiles/nunezlab/D_suzukii_resources/reads/KY_NTeets_2020_2021
fa=/netfiles/nunezlab/D_suzukii_resources/genomes/LBDM_Dsuz_2.1.pri/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.fna

### RUN BWA MEM 2
$bwa mem \
-t $CPU \
$fa \
$rloc/KY${SAMP}_R1.fastq.gz \
$rloc/KY${SAMP}_R2.fastq.gz \
> KY${SAMP}.sam

###
mkdir flagstats
samtools flagstat \
--threads $CPU \
KY${SAMP}.sam \
> flagstats/KY${SAMP}.txt

####
samtools view \
-b \
-q $QUAL \
--threads $CPU  \
KY${SAMP}.sam \
> KY${SAMP}.bam

####
java -Xmx$JAVAMEM \
-jar $PICARD SortSam \
I=KY${SAMP}.bam \
O=KY${SAMP}.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

####
mkdir dupstats
java -Xmx$JAVAMEM \
-jar $PICARD MarkDuplicates \
I=KY${SAMP}.srt.bam \
O=KY${SAMP}.srt.rmdp.bam \
M=dupstats/KY${SAMP}.dupstat.txt \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

####
mkdir KY${SAMP}_qual

$qualimap bamqc \
-bam KY${SAMP}.srt.rmdp.bam  \
-outdir  KY${SAMP}_qual \
--java-mem-size=$JAVAMEM

echo "done"