#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=consensus_pt1

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=28:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

# Submit job array
#SBATCH --array=2-560

# Name output of this job using %x=job-name and %j=job-id
#SBATCH -o ./slurmOutput/consensus_pt1.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------

# This script will run the first half of the consensus script. 

#--------------------------------------------------------------------------------

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly
#WORKING_FOLDER_NETFILES=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

#--------------------------------------------------------------------------------

### import master partition file
guide=$WORKING_FOLDER_SCRATCH/consensus/dat.win.partitions.txt

echo ${SLURM_ARRAY_TASK_ID}

init_bck=$(cat ${guide} | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
final_bck=$(cat ${guide} | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo $init_bck $final_bck
shiftn=`expr $init_bck - 1`
echo $shiftn

#var1=$(echo "^>Backbone_${init_bck}$" ) 
#var2=$(echo "^>Backbone_${final_bck}$" ) 
#echo $var1 $var2

###please run as
#-># sbatch --reservation=jcnunez_test  run_consesus_caller.sh
###please run as

hostname

# Show the ulimit at startup
echo "ulimit at startup"
ulimit -n
# Change the limit
#ulimit -n 1048576
# Show new value
#echo "ulimit after change"
#ulimit -n

# change some stuff
chmod 777 *

# load modules
module load blasr
module load python/python-2.7.5

#(1) backbone_raw.fasta by DBG2OLC
#backbone=/gpfs2/scratch/jcnunez/barnacle_genome/DBG2LC/KCov2_Min30_Adp0.02/KCov2_Min30_Adp0.02.backbone_raw.fasta
#backbone=Backbone_1.fa
#(2) DBG2OLC_Consensus_info.txt by DBG2OLC
#cons_info=/gpfs2/scratch/jcnunez/barnacle_genome/DBG2LC/KCov2_Min30_Adp0.02/DBG2OLC_Consensus_info.txt
#cons_info=test.txt

###
#(3) DBG contigs (in fasta format)
#(4) PacBio reads (in fasta format)
#ONT_FA=/gpfs2/scratch/jcnunez/barnacle_genome/DBG2LC/ONTP2_phase1.pass.FQtoFA.fasta
#Contigs=/gpfs2/scratch/jcnunez/barnacle_genome/sparse_Assb/array_contigs/N2_E1.Contigs.txt
#cat $Contigs $ONT_FA > ctg_pb.fasta

####
mkdir consensus_dir_chunked_Jan9

### Run consensus
sh split_and_run_sparc.pt1.sh \
gen_chunks/gen_chunks.$init_bck.$final_bck.fasta \
chunks/chunk.$init_bck.$final_bck.txt \
ctg_pb.fasta \
./consensus_dir_chunked_Jan9 \
2 \
32 \
$shiftn > cns_log.txt 2>&1

echo "done"