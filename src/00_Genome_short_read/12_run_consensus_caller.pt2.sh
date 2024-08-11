#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=consensus_pt2

# Specify partition
#SBATCH --partition=bigmem

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=200G

# Request CPUs
#SBATCH --cpus-per-task=4

# Submit job array
#SBATCH --array=68-69 #1-559%15

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/consensus_pt2.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will run the second half of the consensus script. 

#--------------------------------------------------------------------------------

# Load modules
module load blasr
spack load python@2.7.18

#--------------------------------------------------------------------------------

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly
WORKING_FOLDER_NETFILES=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

#--------------------------------------------------------------------------------
# Change to consensus directory
#cd $WORKING_FOLDER_SCRATCH/consensus

# Input files for consensus: 
#(1) backbone_raw.fasta by DBG2OLC
#backbone=$WORKING_FOLDER_SCRATCH/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_KmC_2_MinOv_100_Adth_0.01.backbone_raw.fasta
#(2) DBG2OLC_Consensus_info.txt by DBG2OLC
#cons_info=$WORKING_FOLDER_SCRATCH/DBG2OLC/DBG2OLC_KmC_2_MinOv_100_Adth_0.01/DBG2OLC_Consensus_info.txt
#(3) DBG contigs (in fasta format)
#Contigs=$WORKING_FOLDER_NETFILES/SparseAssembler/SparseAssembler_101_2_1/Contigs.txt
#(4) ONT reads (in fasta format) - converted in step 9 part 1
#ONT_FA=$WORKING_FOLDER_SCRATCH/consensus/Nuc.2000.fltlong.FQtoFA.fasta

# Cat contigs and the raw reads for consensus - cat in step 9 part 2
#cat $Contigs $ONT_FA > ctg_ont.fasta

#--------------------------------------------------------------------------------

### import master partition file (first row is column headers)
guide=$WORKING_FOLDER_SCRATCH/consensus/dat.win.partitions.txt

echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, delete the first row (i.e., the column headers), print the first column, then extract the row based on the Slurm array task ID
init_bck=$(cat ${guide} | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
final_bck=$(cat ${guide} | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo $init_bck $final_bck
shiftn=`expr $init_bck - 1`
echo $shiftn

#--------------------------------------------------------------------------------
# Move to working directory
cd $WORKING_FOLDER_SCRATCH/consensus

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "consensus_dir_chunked_Aug2024" ]
then echo "Working consensus_dir_chunked_Aug2024 folder exist"; echo "Let's move on."; date
else echo "Working consensus_dir_chunked_Aug2024 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/consensus/consensus_dir_chunked_Aug2024; date
fi

if [ -d "final_assembly" ]
then echo "Working final_assembly folder exist"; echo "Let's move on."; date
else echo "Working final_assembly folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/consensus/final_assembly; date
fi

#--------------------------------------------------------------------------------

# We need to open a lot of files to distribute the above file into lots of smaller files

# Show the ulimit at startup
#echo "ulimit at startup"
#ulimit -n
# Change the limit
#ulimit -n 1048576
# Show new value
#echo "ulimit after change"
#ulimit -n

# change permissions for consensus
chmod 777 *


# change permissions for accompanying scripts
cd /gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra
chmod 777 *

#--------------------------------------------------------------------------------

# Move to working directory
cd $WORKING_FOLDER_SCRATCH/consensus

# Run consensus

### Run consensus (i.e. run split_and_run_sparc.pt2.sh)
sprun_pt2=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra/split_and_run_sparc.pt2.testrun.sh

$sprun_pt2 \
gen_chunks/gen_chunks.${init_bck}.${final_bck}.fasta \
chunks/chunk.${init_bck}.${final_bck}.txt \
ctg_ont.fasta \
$WORKING_FOLDER_SCRATCH/consensus/consensus_dir_chunked_Aug2024 \
2 \
32 \
$shiftn > cns_log_pt2.txt 2>&1
# 2>&1 redirects stderr to stdout

echo "done"