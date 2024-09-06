#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=busco_ntlink

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Request CPUs per task
#SBATCH -c 1

# Reserve walltime -- hh:mm:ss --30 hour limit
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/busco.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name pilon #create and name the environment
conda activate pilon #activate the environment
#conda install bioconda::pilon # install the program

#--------------------------------------------------------------------------------











#--------------------------------------------------------------------------------

conda deactivate