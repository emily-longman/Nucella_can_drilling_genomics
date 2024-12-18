#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SNPeff

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# SNPeff code

# Create directory 
#cd /netfiles/pespenilab_share/Nucella/processed
#mkdir N.can_genome_Dec2024

# Re-save genome as sequences.fa then gzip
#awk '/^>/ {$0=$1} 1' N.canaliculata_assembly.fasta.softmasked > sequences.fa
#gzip sequences.fa

# Move gtf file, rename and gzip
#scp braker.gtf /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/
#mv braker.gtf genes.gtf

# Use gtf to create gff 
#perl /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/24_gtf_to_gff.pl \
#< /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/braker/braker_cDNA/braker/braker.gtf \
#-o myfile.gff

# Move gff and rename
#scp myfile.gff /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/
#mv myfile.gff N.can.gff


# Create protein file
#module load singularity
#cd /netfiles/nunezlab/Shared_Resources/Software/AGAT
#singularity run agat_1.0.0--pl5321hdfd78af_0.sif
#agat_sp_extract_sequences.pl --gff /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.can.gff \
#-f /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.canaliculata_assembly.fasta.softmasked -p -o /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/N.can.protein.fa

# gzip protein file
#gzip N.can.protein.fa

#--------------------------------------------------------------------------------

# Build directory

# Load SNPeff
module load snpeff

# Move SNPeff config file over
datdir=/netfiles/pespenilab_share/Nucella/processed
param=/netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/snpEff.config

snpeff build -dataDir $datdir -c $param  -gtf22 -v N.can_genome_Dec2024

#--------------------------------------------------------------------------------

# Running 

#module load snpeff

#vcf=...
#datdir=/netfiles/pespenilab_share/Nucella/processed
#param=/netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/snpEff.config

#snpeff -c $param -dataDir $datdir N.can_genome_Dec2024 $vcf > file.ann.vcf