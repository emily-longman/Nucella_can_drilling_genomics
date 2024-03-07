#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=fast_QC 

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --ntasks-per-node=5 # this is number of CPUs you want to use for parallel computing [also referred to as threads] - note not all programs will allow for parallelism, but if they do then its good to use as it helps your jobs run faster

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/fastqc.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

set -e
set -x #useful for debugging your script. It causes each command in the file to be printed to the log file as it is executed, with a + in front of it

# Call fastp package (you would change this to just say 'module load [fastQC version], as I imagine you won't need to install with conda )
module load fastqc-0.11.7-gcc-7.3.0-vcaesw7

#Set up directory
RUN_PATH="/netfiles/pespenilab_share/Nucella/raw/Shortreads/082123-XP-fc1-L7-I8I8-18048-181-396223829/082123-XP-fc1-L7-I8I8-18048-181-688416898/082123-XP-fc1-ds.0852b8f70cc6420995dc0100cb059460/082123-XP-fc1-L7-18048-181-08232023"
cd $RUN_PATH

# create loop to run fastQC on each input file

for file in $(ls -1 _*.gz) #this is listing all your raw data files (assuming they end in .gz)
do
    SAMPLE=`basename $file` #keep the same prefix naming 
    fastqc -t 5 ${SAMPLE} -o /gpfs1/home/e/l/elongman/scratch/Nucella_can_drilling_genomics/data/processed/fastqc #run fastQC, -t here refers to threads, and as we said 5 tasks per node above, we say 5 here, -o indicates output directory
done

