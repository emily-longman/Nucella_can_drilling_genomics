#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=scp.reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --ntasks-per-node=1 # this is number of CPUs you want to use for parallel computing [also referred to as threads] 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=5:00:00 #<= this may depend on your resources

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G #<= this may depend on your resources

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email upd#ates to

#--------------------------------------------------------------------------------

# Move reads from multiple runs into one folder 

scp /netfiles/pespenilab_share/Nucella/raw/Shortreads/082123-XP-fc1-L7-I8I8-18048-181-396223829/082123-XP-fc1-L7-I8I8-18048-181-688416898/082123-XP-fc1-ds.0852b8f70cc6420995dc0100cb059460/082123-XP-fc1-L7-18048-181-08232023/*fastq.gz /netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads/
scp /netfiles/pespenilab_share/Nucella/raw/Shortreads/082123-XP-fc1-L8-I8I8-18048-181-396225832/082123-XP-fc1-L8-I8I8-18048-181-688550906/082123-XP-fc1-ds.666cbe3e3df14e29a6313b5959b79319/082123-XP-fc1-L8-18048-181-08232023/*fastq.gz /netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads/  
scp /netfiles/pespenilab_share/Nucella/raw/Shortreads/Test_run/071023-x-fc1-L2-I8I8-18048-172-07172023-HS-393458067/071023-x-fc1-L2-I8I8-18048-172-07172023-682298951/071023-x-fc1-ds.1dacbba9ea27424790a923cf94e69163/071023-x-fc1-L2-18048-172-07172023/*fastq.gz /netfiles/pespenilab_share/Nucella/raw/Shortreads/All_shortreads/ 

