#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=SFS_sites

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=28:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=80G 

# Submit job array
#SBATCH --array=0-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/Diversity_sites.%A_%a.out # Standard output

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

#Name of pipeline
PIPELINE=Diversity_stats_sites

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Begin Pipeline

# This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "Your unique run id is:" $unique_run_id

echo $PIPELINE
echo $WORKING_FOLDER

if [[ -e "${PIPELINE}.warnings.log" ]]
then echo "Warning log exist"; echo "Let's move on."; date
else echo "Warning log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/${PIPELINE}.warnings.log; date
fi

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "site_frequency_spectrum" ]
then echo "Working site_frequency_spectrum folder exist"; echo "Let's move on."; date
else echo "Working site_frequency_spectrum folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/site_frequency_spectrum; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/site_frequency_spectrum

#--------------------------------------------------------------------------------

# Establish array

arr=("FB" "HC" "MP")
L="${arr[$SLURM_ARRAY_TASK_ID]}"
echo $L

#--------------------------------------------------------------------------------

# Estimating the Site Frequency Spectrum (SFS)

#Estimation of the SFS for all sites using the FOLDED SFS
realSFS ${INPUT}/${L}_GL.saf.idx \
-maxIter 1000 \
-tole 1e-6 \
-P 1 \
-fold 1 \
> ${OUTPUT}/${L}.sfs

#--------------------------------------------------------------------------------

# Estimate theta diversity stats

# Estimate the thetas for each site
realSFS saf2theta ${INPUT}/${L}_GL.saf.idx \
-sfs ${OUTPUT}/${L}.sfs \
-outname ${OUTPUT}/${L}

# Estimate thetas using the SFS
thetaStat do_stat ${OUTPUT}/${L}.thetas.idx

# Estimate thetas using the SFS on a sliding window
thetaStat do_stat ${OUTPUT}/${L}.thetas.idx \
-win 50000 \
-step 10000 \
-outnames ${OUTPUT}/${L}.thetasWindow.gz

# Cut the first column becuase formatted a bit funny
cut -f2- ${OUTPUT}/${L}.thetas.idx.pestPG > ${OUTPUT}/${L}.thetas

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${L} " completed" >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "pipeline completed" $(date)