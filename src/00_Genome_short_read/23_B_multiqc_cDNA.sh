#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=multiqc_CDNA

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 # on one node
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=24:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------
 
# Call package (installed with conda)
module load python3.11-anaconda/2023.09-0
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name multiqc #If you haven't already done so, create and name the environment
source activate multiqc #activate the environment
#conda install -c bioconda multiqc # f you haven't already done so, install the program
conda activate multiqc 

#--------------------------------------------------------------------------------

#Define important file locations

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly

#RAW cDNA indicates the folder where the raw reads are stored.
RAW_READS=/netfiles/pespenilab_share/Nucella/raw/cDNA/all_cDNA

#--------------------------------------------------------------------------------

cd $WORKING_FOLDER_SCRATCH

# Make Quast directory 
if [ -d "cDNA_multiqc" ]
then echo "Working cDNA_multiqc folder exist"; echo "Let's move on."; date
else echo "Working cDNA_multiqc folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/cDNA_multiqc; date
fi

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER_SCRATCH

# Lets do some QC on the reads
# Run multiqc on all of the reads
multiqc $WORKING_FOLDER_SCRATCH/cDNA_fastqc \
-n multiqc_report_cDNA_all.html \
-o multiqc