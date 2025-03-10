#!/usr/bin/env bash  
#  
#SBATCH -J Index_ref  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x_%j.out    
#SBATCH -p general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Load software 
module load gcc/13.3.0-xp3epyt
module load samtools/1.19.2-pfmpoam
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

#--------------------------------------------------------------------------------

# Set folders and file locations
ref=/netfiles/thermofly/GENOMES/basisetae/GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa

#--------------------------------------------------------------------------------

cd /netfiles/thermofly/GENOMES/basisetae

# Index reference

# bwa mem
$bwa index $ref 

# picard
java -jar $PICARD CreateSequenceDictionary \
R=$ref

# samtools
samtools faidx $ref