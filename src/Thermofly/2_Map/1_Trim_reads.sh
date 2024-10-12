#!/usr/bin/env bash  
#  
#SBATCH -J Trim_reads  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/myarray.%A_%a.out  
#SBATCH -p bluemoon  
#SBATCH --array=1-22

#--------------------------------------------------------------------------------

# Trim reads using fastp

# Load software  
fastp=/netfiles/thermofly/shared_software/fastp

#--------------------------------------------------------------------------------

echo ${SLURM_ARRAY_TASK_ID}

# Set folders and file locations
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
meta=$working_folder/METADATA/Thermofly_metadata.tsv

#--------------------------------------------------------------------------------

# Use metadata file to extract sample names and forward and reverse reads
FIL1=$(cat ${meta} | awk -F '\t' '{print $24}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
FIL2=$(cat ${meta} | awk -F '\t' '{print $25}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

SAMP_NAME=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo $FIL1
echo $FIL2
echo $SAMP_NAME

#--------------------------------------------------------------------------------

# Create output folders
cd $working_folder
mkdir fastp_reports
mkdir cleaned_reads

#--------------------------------------------------------------------------------
# This part of the pipeline will generate log files to record warnings and completion status

echo ${SAMP_NAME} "Trimming reads"

# Call fastp and do some light trimming
$fastp \
-i ${FIL1} \
-I ${FIL2} \
-o $working_folder/cleaned_reads/${SAMP_NAME}_R1_clean.fq.gz \
-O $working_folder/cleaned_reads/${SAMP_NAME}_R2_clean.fq.gz \
--detect_adapter_for_pe \
--trim_front1 12 \
--trim_poly_g \
--thread 6 \
--cut_right \
--cut_right_window_size 6 \
--qualified_quality_phred 20 \
--html $working_folder/fastp_reports/${SAMP_NAME}_clean.html \
--json $working_folder/fastp_reports/${SAMP_NAME}_clean.json

# Adapter trimming is enabled by default using overlap analysis, but for PE data you can also specify adapter sequence auto-detection by specifying --detect_adapter_for_pe
# For read 1 of PE data, the front trimming settings are --trim_front1
# Detect and trim polyG in read tails using --trim_poly_g
# Per read cutting by quality options:
### --cut_right: Move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, then stop
### --cut_right_window_size: The window size option for cut_right
### --qualified_quality_phred: The quality vlue that a base is qualified 

#--------------------------------------------------------------------------------

# Inform completion of pipeline
echo "pipeline completed" $(date)