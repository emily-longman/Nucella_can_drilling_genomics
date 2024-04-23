#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ntLinkNu.10000 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --cpus-per-task=40 
#SBATCH --nodes=1 # on one node

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=28:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=95G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/ntLinkNu.10000.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

# ---------------------------

# Move to the directory where the output files will be saved
cd /netfiles/pespenilab_share/Nucella/processed/Base_Genome/FL_10000

# ---------------------------

# Instead try to just install ntlink
# curl -L --output ntLink-1.3.9.tar.gz https://github.com/bcgsc/ntLink/releases/download/v1.3.9/ntLink-1.3.9.tar.gz && tar xvzf ntLink-1.3.9.tar.gz 
# This didn't work either

# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name ntlink python=3.11.0 #create and name the environment
source activate ntlink #activate the environment
#conda install -c bioconda ntlink # install the program
conda install -c bioconda -c conda-forge ntlink
#conda install -c bioconda --file requirements.txt


link_pool=/netfiles/pespenilab_share/Nucella/raw/ONT/FC_all.ONT.nuc.fastq.gz
asm=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Shasta_assemblies/ShastaRun10000/Assembly.fasta 

# Move the assembly to the Base_genome file before scafolding
cp $asm /netfiles/pespenilab_share/Nucella/processed/Base_Genome/FL_10000/ 

#Command from online#
ntLink scaffold target=Assembly.fasta reads=$link_pool k=32 w=250

#Code from Joaquin
#ntLink_rounds run_rounds_gaps \
#target=Assembly.fasta \
#reads=$link_pool k=32 w=100 t=40 rounds=10

conda deactivate

echo "done"