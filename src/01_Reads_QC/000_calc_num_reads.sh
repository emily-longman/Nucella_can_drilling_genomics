#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Cacl_num_reads

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=5 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=70G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to


#--------------------------------------------------------------------------------

# Calculate the number of reads in each of my HiSeq and NovaSeq X Plus lanes

# Short reads directory
SHORT_READS_DIR=/netfiles/pespenilab_share/Nucella/raw/Shortreads

Test_run=$SHORT_READS_DIR/Test_run/071023-x-fc1-L2-I8I8-18048-172-07172023-HS-393458067/071023-x-fc1-L2-I8I8-18048-172-07172023-682298951/071023-x-fc1-ds.1dacbba9ea27424790a923cf94e69163/071023-x-fc1-L2-18048-172-07172023

#--------------------------------------------------------------------------------

cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "General_stats" ]
then echo "Working General_stats folder exist"; echo "Let's move on."; date
else echo "Working General_stats folder doesnt exist. Let's fix that."; mkdir General_stats; date
fi

#--------------------------------------------------------------------------------

Output_dir=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/General_stats

cd $Test_run

# Calculate the number of reads for R1 in Test run

{
for i in *R1*.gz; do
echo $(cat $i|wc -l)/4|bc
done 
} > $Output_dir/Test_run_summary.R1  

{
for i in *R2*.gz; do
echo $(cat $i|wc -l)/4|bc
done 
} > $Output_dir/Test_run_summary.R2  
