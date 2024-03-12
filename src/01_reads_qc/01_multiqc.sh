#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=multiqc

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --ntasks-per-node=1 # this is number of CPUs you want to use for parallel computing [also referred to as threads] - note not all programs will allow for parallelism, but if they do then its good to use as it helps your jobs run faster

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=2:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#Set up directory (make sure this is where the fastQC outputs are!)
cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastqc


# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name multiqc #create and name the environment
source activate multiqc #activate the environment
conda install -c bioconda multiqc # install the program
conda activate multiqc 

#run multiqc in the directory where all of the individual fastqc outputs are
multiqc .