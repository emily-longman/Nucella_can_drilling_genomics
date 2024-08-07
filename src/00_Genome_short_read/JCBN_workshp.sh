# Load modules
module load blasr
spack load python@2.7.18

#--------------------------------------------------------------------------------

# Working folder is core folder where this pipeline is being run.
WORKING_FOLDER_SCRATCH=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/short_read_assembly
WORKING_FOLDER_NETFILES=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/short_read_assembly

guide=$WORKING_FOLDER_SCRATCH/consensus/dat.win.partitions.txt

echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, delete the first row (i.e., the column headers), print the first column, then extract the row based on the Slurm array task ID
init_bck=$(cat ${guide} | sed '1d' | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
final_bck=$(cat ${guide} | sed '1d' | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo $init_bck $final_bck


shiftn=`expr $init_bck - 1`
echo $shiftn


# Move to working directory
cd $WORKING_FOLDER_SCRATCH/consensus

# Generate Folders and files

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "consensus_dir_chunked_July2024" ]
then echo "Working consensus_dir_chunked_July2024 folder exist"; echo "Let's move on."; date
else echo "Working consensus_dir_chunked_July2024 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER_SCRATCH/consensus/consensus_dir_chunked_July2024; date
fi

sprun=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra/split_and_run_sparc.pt1.sh

$sprun \
gen_chunks/gen_chunks.${init_bck}.${final_bck}.fasta \
chunks/chunk.${init_bck}.${final_bck}.txt \
ctg_ont.fasta \
$WORKING_FOLDER_SCRATCH/consensus/consensus_dir_chunked_July2024 \
2 \
32 \
$shiftn > cns_log.txt 2>&1
# 2>&1 redirects stderr to stdout

echo "done"



###

### Run consensus

spspa2=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/src/00_Genome_short_read/12_consensus_scripts_extra/split_and_run_sparc.pt2.sh

sh $spspa2 \
gen_chunks/gen_chunks.${init_bck}.${final_bck}.fasta \
chunks/chunk.${init_bck}.${final_bck}.txt \
ctg_ont.fasta \
$WORKING_FOLDER_SCRATCH/consensus/consensus_dir_chunked_July2024 \
2 \
32 \
$shiftn > cns_log_step2.txt 2>&1

echo "done"



###

 



cmd="blasr -nproc $ncpus consensus_dir_chunked_July2024/backbone-20841.reads.fasta consensus_dir_chunked_July2024/backbone-20841.fasta  -bestn 1 -m 5 -minMatch 19 -out TEST.mapped.m5"
echo $cmd ; eval $cmd


cmd="Sparc m TEST.mapped.m5 b consensus_dir_chunked_July2024/backbone-20841.fasta k 1 c 2 g 1 HQ_Prefix Contig boost 5 t 0.2 o TEST"
echo $cmd ; eval $cmd     

d=TEST
cmd="mv ${d}.consensus.fasta ${d}.fasta"
echo $cmd ; eval $cmd
