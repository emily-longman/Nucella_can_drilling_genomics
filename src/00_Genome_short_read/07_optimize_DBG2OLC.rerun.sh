#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=DBG2OLC_optimized_rerun

# Specify partition
#SBATCH --partition=bigmemwk

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss --7 day limit
#SBATCH --time=07-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# DBG2OLC executable
DBG2OLC=/gpfs1/home/e/l/elongman/software/DBG2OLC
# Quast executable
quast=/netfiles/nunezlab/Shared_Resources/Software/quast-5.2.0/quast.py

#--------------------------------------------------------------------------------

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

# Long reads
#ONT_FOLDER=/netfiles/pespenilab_share/Nucella/raw/ONT
#If you haven't done it yet, unzip the files 
#gunzip $ONT_FOLDER/FC_all.ONT.nuc.fastq.gz
ONT=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/ONT_fltlong/Nuc.2000.fltlong.fastq

#--------------------------------------------------------------------------------

# Determine parameter combination to process
AT=0.01
MO=100
kCT=2

label=DBG2OLC_KmC_${kCT}_MinOv_${MO}_Adth_${AT}

echo ${label} ${kCT} ${MO} ${AT}

#--------------------------------------------------------------------------------

cd $WORKING_FOLDER/DBG2OLC

# Make directory for each DBG2OLC parameter combination
if [ -d "${label}" ]
then echo "Working ${label} folder exist"; echo "Let's move on."; date
else echo "Working ${label} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/DBG2OLC/${label}; date
fi

# Make Quast directory for each parameter combination
cd $WORKING_FOLDER/DBG2OLC/Quast
if [ -d "${label}" ]
then echo "Working ${label} folder exist"; echo "Let's move on."; date
else echo "Working ${label} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/DBG2OLC/Quast/${label}; date
fi

#--------------------------------------------------------------------------------

# Move to the directory where the output files will be saved
cd $WORKING_FOLDER/DBG2OLC/${label}

# Using the SparseAssembler with the largest N50 (SparseAssembler_101_2_1)
# k= 101, NodeCovTh = 2, EdgeCovTh =1 (the last parameter didn't have an effect)

# Use DBG2OLC to construct short but accurate contigs  
$DBG2OLC \
k 17 \
AdaptiveTh ${AT} \
KmerCovTh ${kCT} \
MinOverlap ${MO} \
RemoveChimera 1 \
Contigs $WORKING_FOLDER/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt \
f $ONT

#--------------------------------------------------------------------------------

# Rename output backbone file with label name

echo "renaming outputs"

mv backbone_raw.fasta ${label}.backbone_raw.fasta

#--------------------------------------------------------------------------------

# Run quast
$quast $WORKING_FOLDER/DBG2OLC/${label}/${label}.backbone_raw.fasta \
-o $WORKING_FOLDER/DBG2OLC/Quast/${label}

#--------------------------------------------------------------------------------

echo "done"