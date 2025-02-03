#!/usr/bin/env bash  
#  
#SBATCH -J SnpEff  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 4:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p general 
#SBATCH --array=0-99%10
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Generate SnpEff annotation

# Load software  
module load snpeff/5.2c

#--------------------------------------------------------------------------------

# Set folders and file locations
genome_folder=/netfiles/thermofly/GENOMES/basisetae
ref=/netfiles/thermofly/GENOMES/basisetae/GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa

#--------------------------------------------------------------------------------

# Create output folders
cd $genome_folder
mkdir D.basisetae_SnpEff

#--------------------------------------------------------------------------------

# Move a copy of the genome, gff, config file (edit file to add D. basisetae at end) 
#scp $ref $genome_folder/D.basisetae_SnpEff
#scp GCA_035041595.1_ASM3504159v1_genomic.fna.out.gff $genome_folder/D.basisetae_SnpEff
#scp /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/snpEff.config $genome_folder/D.basisetae_SnpEff

# Try renaming the genome as 
cd $genome_folder/D.basisetae_SnpEff
scp $ref $genome_folder/D.basisetae_SnpEff
mv GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa D.basisetae_SnpEff.genome

#--------------------------------------------------------------------------------

cd $genome_folder

# Define important file locations
datdir=$genome_folder
param=$genome_folder/D.basisetae_SnpEff/snpEff.config

snpeff build -dataDir $datdir -c $param -gff3 -noCheckCds -v D.basisetae_SnpEff