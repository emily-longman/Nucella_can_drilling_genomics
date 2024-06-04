#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=multiqc_qualimap

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=20G 

# Submit job array
#SBATCH --array=1-576%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will calculate average coverage 

#--------------------------------------------------------------------------------
# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
conda create --name multiqc #create and name the environment
source activate multiqc #activate the environment
conda install -c bioconda multiqc # install the program
conda activate multiqc 
#--------------------------------------------------------------------------------

## PREPARE GUIDE FILES
## Read guide files
# This is a file with the name all the samples to be processed. One sample name per line with all the info.

SAMPLE_FILE=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF/GuideFileMerged.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
## Snail_ID  Sample#   Merged_name 1    Merged_name 2    Merged_name 3   Merged_bam_name
##  FB1-1     S84     FB1-1_S84_L002   FB1-1_S84_L007   FB1-1_S84_L008     FB1-1_S84
##  FB1-2     S173    FB1-2_S173_L002  FB1-2_S173_L007  FB1-2_S173_L008    FB1-2_S173
##  FB1-5     S109    FB1-5_S109_L002  FB1-5_S109_L007  FB1-5_S109_L008    FB1-5_S109
##  ...
##  MP9-9     S191    MP9-9_S191_L002  MP9-9_S191_L007  MP9-9_S191_L008    MP9-9_S191
##  MP9-10    S26     MP9-10_S26_L002  MP9-10_S26_L007  MP9-10_S26_L008    MP9-10_S26

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $6}' $SAMPLE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo $i

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_VCF

cd $WORKING_FOLDER

if [ -d "multiqc_qualimap" ]
then
echo "Working multiqc_qualimap folder exist"
echo "lets move on"
else 
echo "Folder doesnt exist. Lets fix that"
mkdir multiqc_qualimap
fi

#--------------------------------------------------------------------------------

# Put a copy of each qualimap file in the multiqc_qualimap folder
scp $WORKING_FOLDER/Merged_Bams_qualimap/Qualimap_LaneMerged_${i}/qualimapReport.html $WORKING_FOLDER/multiqc_qualimap/qualimapReport_${i}.html

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER/multiqc_qualimap

# Run multiqc in the directory where all of the individual fastqc outputs are
multiqc . -n multiqc_report_bams.html -o .

conda deactivate