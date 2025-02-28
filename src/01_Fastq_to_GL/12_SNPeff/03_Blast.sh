#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Blast

# Specify partition
#SBATCH --partition=general

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

# This script will run SNPeff on a vcf

#Load modules
ncbi=/gpfs1/home/e/l/elongman/software/ncbi-blast-2.16.0+/bin

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "Gene_ontology" ]
then echo "Working Gene_ontology folder exist"; echo "Let's move on."; date
else echo "Working Gene_ontology folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Gene_ontology; date
fi

cd $WORKING_FOLDER/Gene_ontology

if [ -d "N_canaliculata" ]
then echo "Working N_canaliculata folder exist"; echo "Let's move on."; date
else echo "Working N_canaliculata folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Gene_ontology/N_canaliculata; date
fi

if [ -d "uniprot" ]
then echo "Working uniprot folder exist"; echo "Let's move on."; date
else echo "Working uniprot folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Gene_ontology/uniprot; date
fi

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/Gene_ontology/N_canaliculata

# Move a copy of the protein file over to this directory
#scp /netfiles/pespenilab_share/Nucella/processed/N.can_genome_Dec2024/protein.fa .

#--------------------------------------------------------------------------------

# Format N. canaliculata database
#$ncbi/makeblastdb -in protein.fa -dbtype prot

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/Gene_ontology/uniprot

# Download and unzip uniprot database in $WORKING_FOLDER/Gene_ontology directory
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz 
gunzip -v uniref90.fasta.gz                                                         

#--------------------------------------------------------------------------------


# This single line using the blastp command below will compare your transcript fasta file
# (-query) to the already formatted uniref90 database (-db).
# You can enter 'blastp --help' for a list of the parameters.
# We choose the tab-delimited output format (6) and to only help the top hit (-max_target_seqs)
# and only if it has a minimum evalue of 0.00001.

#blastp -query /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds \
#       -db /data/archive/databases/uniref90/uniprot_uniref90.trinotate.pep \
#       -out /data/project_data/assembly/blast/blastp_vs_uniref90.outfmt6 \
#       -outfmt 6 \
#       -evalue 1e-5 \
#       -max_target_seqs 1