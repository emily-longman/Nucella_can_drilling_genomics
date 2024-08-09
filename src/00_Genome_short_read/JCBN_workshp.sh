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



###########################

#Emily test

backbone_fasta=gen_chunks/gen_chunks.${init_bck}.${final_bck}.fasta 
consensus_fasta=chunks/chunk.${init_bck}.${final_bck}.txt
reads_fasta=ctg_ont.fasta
split_dir=$WORKING_FOLDER_SCRATCH/consensus/consensus_dir_chunked_Aug2024
iterations=2
ncpus=32
nshift=$shiftn 
#> cns_log_pt2.txt 2>&1

echo $backbone_fasta
echo $consensus_fasta
echo $reads_fasta
echo $split_dir
echo $iterations
echo $ncpus
echo $nshift

chunk=backbone-10004

cmd=""

d="${split_dir}/backbone-10004"
#d="${split_dir}/${chunk}"
cmd="blasr -nproc $ncpus ${d}.reads.fasta ${d}.fasta -bestn 1 -m 5 -minMatch 19 -out ${chunk}.mapped.m5"
echo $cmd ; eval $cmd ;
cmd="Sparc m ${chunk}.mapped.m5 b ${d}.fasta k 1 c 2 g 1 HQ_Prefix Contig boost 5 t 0.2 o ${chunk}"
echo $cmd ; eval $cmd ;

# rename and move to final assembly directory
cmd="mv ${chunk}.consensus.fasta $WORKING_FOLDER_SCRATCH/consensus/final_assembly/${chunk}.consensus.fasta"
echo $cmd ; eval $cmd ;

cmd="rm ${chunk}.mapped.m5"
echo $cmd
eval $cmd
cmd="rm ${split_dir}/${chunk}.reads.fasta"
echo $cmd
eval $cmd





for file in $(find ${split_dir} -name "*.reads.fasta"); do
chunk=`basename $file .reads.fasta`

cmd=""
for iter in `seq 1 ${iterations}`; do
# create this to shorten the subsequent commands to more readable length
# I think I got them all substituted correctly.
d="${split_dir}/${chunk}"
cmd="blasr -nproc $ncpus ${d}.reads.fasta ${d}.fasta -bestn 1 -m 5 -minMatch 19 -out ${d}.mapped.m5"
echo $cmd ; eval $cmd
cmd="Sparc m ${d}.mapped.m5 b ${d}.fasta k 1 c 2 g 1 HQ_Prefix Contig boost 5 t 0.2 o ${d}"
echo $cmd ; eval $cmd        
if [[ ${iter} -lt ${iterations} ]]
then
# rename and move to final assembly directory
cmd="mv ${d}.consensus.fasta $WORKING_FOLDER_SCRATCH/consensus/final_assembly/${d}.consensus.fasta"
echo $cmd ; eval $cmd
fi
done
