#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Blast_blastp

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=7-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run blast 

#Load modules
ncbi=/gpfs1/home/e/l/elongman/software/ncbi-blast-2.16.0+/bin

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics

PROTEIN_FILE=/netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/protein.fa

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/Gene_ontology

#$ncbi/blastp -query test.txt \
#-db $WORKING_FOLDER/Gene_ontology/uniprot/uniref90 \
#-out $WORKING_FOLDER/Gene_ontology/blastp_vs_uniref90.outfmt6 \
#-outfmt 6 \
#-evalue 1e-5 \
#-max_target_seqs 1

# Use the blastp command to compare the Nucella protein file with the uniprot database

$ncbi/blastp -query $PROTEIN_FILE \
-db $WORKING_FOLDER/Gene_ontology/uniprot/uniref90 \
-out $WORKING_FOLDER/Gene_ontology/blastp_vs_uniref90.outfmt6 \
-outfmt 6 \
-evalue 1e-5 \
-max_target_seqs 1

# query <File_In>: Input file name
# db: BLAST database name
# outfmt 6: alignment view options: Tabular
# evalue: Expectation value (E) threshold for saving hits. Default = 10 
# max_target_seqs: Maximum number of aligned sequences to keep 

#blastp -query /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds \
#       -db /data/archive/databases/uniref90/uniprot_uniref90.trinotate.pep \
#       -out /data/project_data/assembly/blast/blastp_vs_uniref90.outfmt6 \
#       -outfmt 6 \
#       -evalue 1e-5 \
#       -max_target_seqs 1