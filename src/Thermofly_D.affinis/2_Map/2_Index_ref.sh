#!/usr/bin/env bash  
#  
#SBATCH -J Index_ref_D.affinis 
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x_%j.out    
#SBATCH -p bluemoon  
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Load software 
spack load samtools@1.10
bwa=/netfiles/nunezlab/Shared_Resources/Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar

#bwa=/netfiles/thermofly/shared_software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.avx2
#PICARD=/netfiles/thermofly/shared_software/picard/build/libs/picard.jar

#--------------------------------------------------------------------------------

# Set folders and file locations
ref=/netfiles/thermofly/GENOMES/supergenome_affinis/supergenome.aff.alg.atha.masked.fasta

#--------------------------------------------------------------------------------

cd /netfiles/thermofly/GENOMES/affinis

# Index reference

# bwa mem
$bwa index $ref 

# picard
java -jar $PICARD CreateSequenceDictionary \
R=$ref

# samtools
samtools faidx $ref