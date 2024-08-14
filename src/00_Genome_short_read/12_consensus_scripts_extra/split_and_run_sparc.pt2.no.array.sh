#!/usr/bin/env bash

###
# USAGE: ./split_and_run_sparc.sh [BACKBONE_FASTA] [CONSENSUS_FASTA] [READS_FASTA] [OUTPUT_DIR] [ITERATIONS] ###

backbone_fasta=$1 #backbone broken into chunks
consensus_fasta=$2 #text file with contig names in chunks
reads_fasta=$3 #contigs and raw ONT reads catted together
split_dir=$4 #file output directory 
iterations=$5 #set to 2
ncpus=$6 #set to 32
nshift=$7 #shiftn

echo $backbone_fasta
echo $consensus_fasta
echo $reads_fasta
echo $split_dir
echo $iterations
echo $ncpus
echo $nshift


# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly
WORKING_FOLDER_NETFILES=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

### need scripts
##split_reads_by_backbone.py

### Make some aliaces
#alias blasr="/gpfs1/home/j/c/jcnunez/software/blasrmc/alignment/bin/blasrmc"
#alias Sparc="/gpfs1/home/e/l/elongman/software/Sparc"
Sparc=/gpfs1/home/e/l/elongman/software/Sparc

#clean the directory first  --- no longer required
#find ${split_dir} -name "backbone-*" -delete

#===#/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra/split_reads_by_backbone.py -b ${backbone_fasta} -o ${split_dir} -r ${reads_fasta} -c ${consensus_fasta} -s ${nshift}

for file in $(find ${split_dir} -name "*.reads.fasta"); do
chunk=`basename $file .reads.fasta`

cmd=""
for iter in `seq 1 ${iterations}`; do

d="${split_dir}/${chunk}"

echo $iter
echo $chunk
echo $d
echo "Align contig with blasr and Sparc."

blasr -nproc $ncpus ${split_dir}/${chunk}.reads.fasta ${split_dir}/${chunk}.fasta -bestn 1 -m 5 -minMatch 19 -out $WORKING_FOLDER_SCRATCH/consensus/${chunk}.mapped.m5   
$Sparc m $WORKING_FOLDER_SCRATCH/consensus/${chunk}.mapped.m5 b ${split_dir}/${chunk}.fasta k 1 c 2 g 1 HQ_Prefix Contig boost 5 t 0.2 o $WORKING_FOLDER_SCRATCH/consensus/${chunk}

if [[ ${iter} -lt ${iterations} ]]
then

# Rename and move to final assembly directory
cmd="mv $WORKING_FOLDER_SCRATCH/consensus/${chunk}.consensus.fasta $WORKING_FOLDER_SCRATCH/consensus/final_assembly_no_array/${chunk}.consensus.fasta"
echo $cmd ; eval $cmd
#fi
#done

echo $cmd
eval $cmd


#to save space
cmd="rm ${chunk}.mapped.m5"
echo $cmd
eval $cmd
#cmd="rm ${split_dir}/${chunk}.reads.fasta"
#echo $cmd
#eval $cmd

done

#for confile in $(find ${split_dir} -name "*.consensus.fasta"); do
#cmd="cat ${confile};"
#eval $cmd
#done > $WORKING_FOLDER_SCRATCH/consensus/final_assembly_chunked/final_assembly.${SLURM_ARRAY_TASK_ID}.fasta