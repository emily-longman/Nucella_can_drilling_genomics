#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SparseAssembler

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name SparseAssembler #create and name the environment
source activate SparseAssembler #activate the environment
conda install -c bioconda sparseassembler # install the program
conda activate SparseAssembler 

#--------------------------------------------------------------------------------

# If you haven't done it yet, gunzip the files 
gunzip /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/fastp/*fastq.gz

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly/SparseAssembler

# Use SparseAssembler to construct short but accurate contigs  
SparseAssembler \
LD 0 k 51 g 15 \
NodeCovTh 1 \
EdgeCovTh 0 \
GS 2500000000 \
i1 NC3_R1_cleaned.fastq \
i2 NC3_R2_cleaned.fastq

# In future try k = 31 and k = 71

#--------------------------------------------------------------------------------

# Inform that SparseAssembler is done

echo "done"

conda deactivate