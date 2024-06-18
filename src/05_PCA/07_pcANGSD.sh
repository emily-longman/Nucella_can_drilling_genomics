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

#Load python version 
module load python3.11-anaconda/2023.09-0

# Extra code notes
#pip install Cython #don't need to do because Cython is installed in this python version
#python setup.py build_ext --inplace #needed to rebuild pcangsd after Csenge made 'build' executable

# Installing on my own in /gpfs1/home/e/l/elongman/software (still run into an error when do this)
#module load python3.11-anaconda/2023.09-0
#git clone https://github.com/Rosemeis/pcangsd.git
#cd pcangsd/
#python setup.py build_ext --inplace
#pip install --user -r requirements.txt
# Make pcangsd executable
#pcangsd=/gpfs1/home/e/l/elongman/software/pcangsd/pcangsd/pcangsd.py

# Make pcangsd executable (Csenge version)
pcangsd=/netfiles/pespenilab_share/pcangsd/pcangsd.py

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#Input folder is genotype likelihoods from ANGSD
INPUT=$WORKING_FOLDER/genotype_likelihoods_SNPs

#Name of pipeline
PIPELINE=pcangsd

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

#--------------------------------------------------------------------------------
# Make a copy of the list of bam files for all the Nucella samples and place in the output directory. You'll need this later for making figures.

cp ${INPUT}/Nucella_bam.list ${OUTPUT}

# Then, run PCA and admixture scores with pcangsd:
SUFFIX="Nucella_poly_covmatrix"

python /netfiles/pespenilab_share/pcangsd/pcangsd.py \
-b ${INPUT}/Nucella_SNPs.beagle.gz \
-o ${OUTPUT}/${SUFFIX} \
-threads $CPU 

# PCAngsd accepts either genotype likelihoods in Beagle format generated from BAM files using ANGSD
# -e : Manually select the number of eigenvalues to use in the modelling of individual allele frequencies (the number of clusters is K-1)
# --admix : Individual admixture proportions and ancestral allele frequencies can be estimated assuming K ancestral populations using an accelerated mini-batch NMF method.

# https://www.popgen.dk/software/index.php/PCAngsd