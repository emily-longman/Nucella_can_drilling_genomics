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

# Generate SnpEff database

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

# Move a copy of the genome, gff, config file to directory
#scp $ref $genome_folder/D.basisetae_SnpEff
#scp $genome_folder/Dbas.genepred.hardmask.gff3 $genome_folder/D.basisetae_SnpEff
#scp /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/snpEff.config $genome_folder/D.basisetae_SnpEff

# Note: edit config file to add D. basisetae directory at end

# Rename the genome as seqeuences.fa
cd $genome_folder/D.basisetae_SnpEff
scp $ref $genome_folder/D.basisetae_SnpEff
mv GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa sequences.fa

# Rename the gff files as genes.gff
mv Dbas.genepred.hardmask.gff3 genes.gff
#--------------------------------------------------------------------------------

# Create protein file using AGAT (RHEL6)
#module load singularity
#cd /netfiles/nunezlab/Shared_Resources/Software/AGAT
#singularity run agat_1.0.0--pl5321hdfd78af_0.sif
#agat_sp_extract_sequences.pl --gff /netfiles/thermofly/GENOMES/basisetae/D.basisetae_SnpEff/genes.gff \
#-f /netfiles/thermofly/GENOMES/basisetae/D.basisetae_SnpEff/sequences.fa -p -o /netfiles/thermofly/GENOMES/basisetae/D.basisetae_SnpEff/protein.fa

#--------------------------------------------------------------------------------

cd $genome_folder

# Define important file locations
datdir=$genome_folder
param=$genome_folder/D.basisetae_SnpEff/snpEff.config

snpeff build -dataDir $datdir -c $param -gff3 -noCheckCds -v D.basisetae_SnpEff
