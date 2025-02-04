#!/usr/bin/env bash  
#  
#SBATCH -J SnpEff  
#SBATCH -c 1  
#SBATCH -N 1 # on one node  
#SBATCH -t 4:00:00   
#SBATCH --mem 40G   
#SBATCH --output=./slurmOutput/%x_%j.out 
#SBATCH -p general 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run SNPeff on a vcf

# Load software  
module load snpeff/5.2c

#--------------------------------------------------------------------------------

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
vcf=/netfiles/thermofly/ANALYSES/Dbas_paper/D.basisetae.GATK.pipe.vcf.gz
datdir=/netfiles/thermofly/GENOMES/basisetae
param=/netfiles/thermofly/GENOMES/basisetae/D.basisetae_SnpEff/snpEff.config

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir SnpEff

#--------------------------------------------------------------------------------

# Change directory
cd $working_folder/SnpEff

# Run SnpEff
snpeff -c $param -dataDir $datdir D.basisetae_SnpEff $vcf > $working_folder/SnpEff/D.bas_SNPs.annotate.vcf