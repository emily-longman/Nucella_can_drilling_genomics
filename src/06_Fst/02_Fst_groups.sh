#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Fst_groups

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=05-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G 

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

#Name of pipeline
PIPELINE=Fst_groups

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Begin Pipeline

# Welcome message
echo "Your unique run id is:" $unique_run_id

echo $PIPELINE
echo $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "Fst_groups" ]
then echo "Working Fst_groups folder exist"; echo "Let's move on."; date
else echo "Working Fst_groups folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/Fst_groups; date
fi

#Output folder
OUTPUT=$WORKING_FOLDER/Fst_groups

#--------------------------------------------------------------------------------

cd ${OUTPUT}

# Estimate Fst between my red spruce pop and black spruce:

realSFS ${INPUT}/Drilled_GL_allsites.saf.idx ${INPUT}/NotDrilled_GL_allsites.saf.idx \
-P $CPU \
-nSites \
>${OUTPUT}/Drilled_NotDrilled.sfs

# If you have .saf file larger than -nSites (you can check the number of sites in the .saf.pos file), then the program will loop over the genome and output the results for each block. So each line in your Whit.saf.ml, is an SFS for a region.

#realSFS fst index \
#${INPUT}/Drilled_GL_allsites.saf.idx \
#${INPUT}/NotDrilled_GL_allsites.saf.idx \
#-sfs ${OUTPUT}/Drilled_NotDrilled.sfs \
#-fstout ${OUTPUT}/Drilled_NotDrilled \
#-whichFst 1

#realSFS fst stats ${OUTPUT}/Drilled_NotDrilled.fst.idx 

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part will notify the completion of the pipeline 

echo "pipeline completed" $(date)