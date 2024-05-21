#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Genotype_likelihoods

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------


#Load modules 
module load angsd-0.933-gcc-7.3.0-4wsdzjw
module load samtools-1.10-gcc-7.3.0-pdbkohx

# Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_May2024/Assembly.fasta.k24.w150.z1000.ntLink.8rounds.fa

#Input folder is genotype likelihoods from ANGSD
INPUT=$WORKING_FOLDER/genotype_likelihoods

#--------------------------------------------------------------------------------
# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "site_frequency_spectrum" ]
then echo "Working site_frequency_spectrum folder exist"; echo "Let's move on."; date
else echo "Working site_frequency_spectrum folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/site_frequency_spectrum; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/site_frequency_spectrum
#--------------------------------------------------------------------------------

# Estimating the Site Frequency Spectrum (SFS)

# File suffix to distinguish analysis choices
SUFFIX="SFS"

#Estimation of the SFS for all sites using the FOLDED SFS
realSFS ${INPUT}/Nucella_GL.saf.idx \
-maxIter 1000 \
-tole 1e-6 \
-P 1 \
> ${OUTPUT}/Nucella_${SUFFIX}.sfs

#--------------------------------------------------------------------------------

# Estimate theta diversity stats

# Estimate the thetas for each site
realSFS saf2theta ${OUTPUT}/Nucella.saf.idx \
-sfs ${OUTPUT}/Nucella.sfs \
-outname ${OUTPUT}/Nucella

# Estimate thetas using the SFS
thetaStat do_stat ${OUTPUT}/Nucella.thetas.idx

# Estimate thetas using the SFS on a sliding window
thetaStat do_stat ${OUTPUT}/Nucella.thetas.idx \
-win 50000 \
-step 10000 \
-outnames ${OUTPUT}/Nucella.thetasWindow.gz

# Cut the first column becuase formatted a bit funny
#cut -f2- ${OUTPUT}/Nucella.thetas.idx.pestPG > ${OUTPUT}/Nucella.thetas