#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=pcangsd

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss
#SBATCH --time=2:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Install pcangsd
#spack find
module load python3.11-anaconda
module load gcc/10.5.0

# Set a directory for your project and the name of your virtual environment:
export project_dir=/gpfs1/home/e/l/elongman/venv
export venv_name=pcangsd-venv

# Change into the project directory, create and activate your virtual environment:
#cd $project_dir
#python -m venv $venv_name  
#source $venv_name/bin/activate

# Download PCAngsd via git clone and change to that directory:
#git clone https://github.com/Rosemeis/pcangsd.git
#cd pcangsd/

# Install PCAngsd and run the help command
#pip3 install .
#pcangsd -h

# Change into the project directory and activate your virtual environment:
cd $project_dir
source $venv_name/bin/activate

#--------------------------------------------------------------------------------

# Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

# Scripts folder
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

# Path to bam list.
BAM_LIST=$WORKING_FOLDER/guide_files/Nucella_bam.list

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#--------------------------------------------------------------------------------

# Prepare variables 

# Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "pcangsd" ]
then echo "Working pcangsd folder exist"; echo "Let's move on."; date
else echo "Working pcangsd folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pcangsd; date
fi


#--------------------------------------------------------------------------------

# Make a copy of the list of bam files for all the Nucella samples and place in the output directory. You'll need this later for making figures.
cp $BAM_LIST $WORKING_FOLDER/pcangsd

# Then, run PCA and admixture scores with pcangsd:
echo "Analyse covariance matrix on all individuals"

pcangsd \
-b $WORKING_FOLDER/genotype_likelihoods_all/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6.beagle.gz \
-o $WORKING_FOLDER/pcangsd/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mindepth"$MIN_DEPTH"_maxdepth"$MAX_DEPTH_FACTOR"_pval1e6 \
-t $CPU 

# PCAngsd accepts either genotype likelihoods in Beagle format generated from BAM files using ANGSD
# -e : Manually select the number of eigenvalues to use in the modelling of individual allele frequencies (the number of clusters is K-1)
# --admix : Individual admixture proportions and ancestral allele frequencies can be estimated assuming K ancestral populations using an accelerated mini-batch NMF method.

echo "Transform covariance matrix into PCA"
COV_MAT=$WORKING_FOLDER/pcangsd/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".cov

#--------------------------------------------------------------------------------

# Deactivate the virtual environment

deactivate