#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=angsd_admix

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10

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

#Load software
NGSadmix=/gpfs1/home/e/l/elongman/software/NGSadmix

#--------------------------------------------------------------------------------

#Define important file locations

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL

#This is the location where the reference genome and all its indexes are stored.
REFERENCE=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Aug2024/backbone_raw.fasta

#Scripts folder
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/01_Fastq_to_GL

#--------------------------------------------------------------------------------

# Define parameters
NB_CPU=10 #change accordingly in SLURM header
echo "using #CPUs ==" $NB_CPU

#Use config file (this means you dont need to directly input minimum individual/depth parameters)
source $SCRIPT_FOLDER/03_Call_SNPs/01_config.sh

#Min number of pop to consider for NGS admix
K_MIN=1
#Maximum number of pop to consider for NGS admix
K_MAX=5

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "ngs_admix" ]
then echo "Working ngs_admix folder exist"; echo "Let's move on."; date
else echo "Working ngs_admix folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/ngs_admix; date
fi

cd $WORKING_FOLDER/ngs_admix

if [ -d "K_output" ]
then echo "Working K_output folder exist"; echo "Let's move on."; date
else echo "Working K_output folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/ngs_admix/K_output; date
fi

#--------------------------------------------------------------------------------

## Run admixture for Kmin to Kmax (these set in the config file) 10 times each
for j in {1..10}; 
do

for i in $(seq $K_MIN $K_MAX)
do 
echo $i
$NGSadmix -P $NB_CPU \
-likes $WORKING_FOLDER/genotype_likelihoods_all/Nucella_SNPs_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth"$MAX_DEPTH_FACTOR".beagle.gz \
-minMaf $MIN_MAF -K $i -o $WORKING_FOLDER/ngs_admix/K_output/Nucella_all_maf_K{$i}_run{$j}
done
done