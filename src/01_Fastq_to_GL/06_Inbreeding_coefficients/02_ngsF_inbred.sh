#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ngsF_inbred

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Request CPU
#SBATCH --cpus-per-task=5

# Submit job array
#SBATCH --array=0-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate inbreeding coefficients with ngsF (https://github.com/fgvieira/ngsF/tree/master)

# Load modules
module load singularity/3.7.1
module load ngstools

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_mask/N.canaliculata_assembly.fasta.masked

# Scripts folder.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=5 #change accordingly in SLURM header
echo "using #CPUs ==" $NB_CPU

#--------------------------------------------------------------------------------

# Establish the array
# This is a file with the names. 
arr=("FB" "HC" "MP")
i="${arr[$SLURM_ARRAY_TASK_ID]}"
echo ${i}

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

# Extract parameters from config file
N_IND=$(wc -l $WORKING_FOLDER/guide_files/${i}_bam.list | cut -d " " -f 1) 
MIN_IND_FLOAT=$(echo "($N_IND * $PERCENT_IND)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*} 
MAX_DEPTH=$(echo "($N_IND * $MAX_DEPTH_FACTOR)" |bc -l)

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "ngsF" ]
then echo "Working ngsF folder exist"; echo "Let's move on."; date
else echo "Working ngsF folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/ngsF; date
fi

# Change directory
cd $WORKING_FOLDER/ngsF

if [ -d "${i}" ]
then echo "Working ${i} folder exist"; echo "Let's move on."; date
else echo "Working ${i} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/ngsF/${i}; date
fi

#--------------------------------------------------------------------------------

# Start pipeline

# Run ngsTools/ngsF for all population listed
# Run in two steps: run a faster preliminary/approximated method then use the output for the slower/main algorithm
# Note: ngsF runs 10 repeats of the analysis to prevent convergence to local maxima

echo "Calculating inbreeding for site/population:" ${i}

# Calculate the number of sites 
NSITES=$((`zcat $WORKING_FOLDER/genotype_likelihoods_by_site/${i}/${i}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_inbreed.mafs.gz | wc -l`-1))
echo $NSITES "sites for collection site/population" ${i} "with" $N_IND "individuals"

# Unzip .glf.gz file if you haven't already done so
#gunzip $WORKING_FOLDER/genotype_likelihoods_by_site/${i}/${i}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_inbreed.glf.gz 

# Preliminary search (i.e., run the approximated method to get some reasonable estimates of the parameters)
singularity run $NGS ngsF \
--n_ind $N_IND --n_sites $NSITES \
--seed 12345 \
--min_epsilon 1e-05 \
--glf $WORKING_FOLDER/genotype_likelihoods_by_site/${i}/${i}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_inbreed.glf \
--out $WORKING_FOLDER/ngsF/${i}/${i}.approx_indF \
--approx_EM --init_values r --n_threads 5
	
# Calc inbreeding (i.e., use the output from the preliminary search as initial values for the main (and slower) algorithm)
singularity run $NGS ngsF \
--n_ind $N_IND --n_sites $NSITES \
--min_epsilon 1e-07 \
--glf $WORKING_FOLDER/genotype_likelihoods_by_site/${i}/${i}_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_inbreed.glf \
--out $WORKING_FOLDER/ngsF/${i}/${i}.indF \
--init_values $WORKING_FOLDER/ngsF/${i}/${i}.approx_indF.pars --n_threads 5