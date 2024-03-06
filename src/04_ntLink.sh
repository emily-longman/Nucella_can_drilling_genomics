#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ntLinkNu 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=28:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=295G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/ntLinkNu.%A_%a.out # Standard output

# Move to the directory where the output files will be saved
cd /netfiles02/pespenilab_share/Nucella/processed/Base_Genome/

# ---------------------------
# NEED TO FINISH UPDATING BELOW - ALSO NEED TO GET miniCONDA

### UPDATE!
source /gpfs1/home/j/c/jcnunez/miniconda3/etc/profile.d/conda.sh

conda activate ntlink

link_pool=/netfiles/pespenilab_share/Nucella/processed/FC_all.ONT.nuc.fastq.gz
asm=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Assembly.fasta

cp $asm ./

ntLink_rounds run_rounds_gaps \
target=Assembly.fasta \
reads=$link_pool k=32 w=100 t=40 rounds=10

echo "done"