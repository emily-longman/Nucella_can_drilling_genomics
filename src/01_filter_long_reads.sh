#!/usr/bin/env bash
#
#SBATCH -J fltlong
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/fltlong.%A_%a.out # Standard output
#SBATCH -p bluemoon
#SBATCH --array=2-4

filtlong=/gpfs1/home/j/c/jcnunez/software/Filtlong/bin/filtlong
input=./FC_all.ONT.nuc.fastq.gz

arr=(1000 2000 3500 5000)
L="${arr[$SLURM_ARRAY_TASK_ID]}"
echo $L

$filtlong \
--min_length $L \
$input | gzip > Nuc.$L.fltlong.fastq.gz