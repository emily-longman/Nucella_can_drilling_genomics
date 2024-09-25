#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=pcangsd_lower_p

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit 
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=65G 

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

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Input folder is genotype likelihoods from ANGSD
INPUT=$WORKING_FOLDER/genotype_likelihoods_all

#--------------------------------------------------------------------------------
# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "pcangsd" ]
then echo "Working pcangsd folder exist"; echo "Let's move on."; date
else echo "Working pcangsd folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pcangsd; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/pcangsd

##NOTE: Ideally, change the output folder to the results folder in the github. I just scp the output over

#--------------------------------------------------------------------------------

# Make a copy of the list of bam files for all the Nucella samples and place in the output directory. You'll need this later for making figures.
cp ${INPUT}/Nucella_bam_lower_p.list ${OUTPUT}

# Then, run PCA and admixture scores with pcangsd:
SUFFIX="Nucella_poly_covmatrix_lower_p"

echo "Analyse covariance matrix on all individuals"

pcangsd -b ${INPUT}/Nucella_SNPs_all_lower_p.beagle.gz \
-o ${OUTPUT}/${SUFFIX} \
-t $CPU 

# PCAngsd accepts either genotype likelihoods in Beagle format generated from BAM files using ANGSD
# -e : Manually select the number of eigenvalues to use in the modelling of individual allele frequencies (the number of clusters is K-1)
# --admix : Individual admixture proportions and ancestral allele frequencies can be estimated assuming K ancestral populations using an accelerated mini-batch NMF method.

echo "Transform covariance matrix into PCA"
COV_MAT=${OUTPUT}/${SUFFIX}.cov

#--------------------------------------------------------------------------------

# Deactivate the virtual environment

deactivate