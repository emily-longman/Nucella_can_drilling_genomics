#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=partition_genome

# Specify partition
#SBATCH --partition=bluemoon

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1  

# Reserve walltime -- hh:mm:ss 
#SBATCH --time=80:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x_%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL # indicates if you want an email when the job starts, ends, or both
#SBATCH --mail-user=emily.longman@uvm.edu # where to email updates to

#--------------------------------------------------------------------------------



#SBATCH -o ./slurmOutput/part_info.%A_%a.out # Standard output
#SBATCH --array=1-1764


#1764
###import files
cons_info=/gpfs2/scratch/jcnunez/barnacle_genome/DBG2LC/KCov2_Min30_Adp0.02/DBG2OLC_Consensus_info.txt
backbone=/gpfs2/scratch/jcnunez/barnacle_genome/DBG2LC/KCov2_Min30_Adp0.02/KCov2_Min30_Adp0.02.backbone_raw.fasta

### import master partition file
guide=./dat.win.partitions.txt

echo ${SLURM_ARRAY_TASK_ID}


init_bck=$(cat ${guide} | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
final_bck=$(cat ${guide} | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo $init_bck $final_bck

var1=$(echo "^>Backbone_${init_bck}$" ) 
var2=$(echo "^>Backbone_${final_bck}$" ) 
echo $var1 $var2

## make genome chunks
mkdir ./gen_chunks
cat $backbone | \
sed -n "/$var1/, /$var2/p" | \
sed '$d'  > ./gen_chunks/gen_chunks.$init_bck.$final_bck.fasta

## make info chunks
mkdir ./chunks
cat $cons_info | \
sed -n "/$var1/, /$var2/p" | \
sed '$d'  > ./chunks/chunk.$init_bck.$final_bck.txt

echo "done"