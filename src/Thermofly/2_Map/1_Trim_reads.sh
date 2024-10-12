#!/usr/bin/env bash  
#  
#SBATCH -J QC  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/myarray.%A_%a.out  
#SBATCH -p bluemoon  
#SBATCH --array=2-24 

#--------------------------------------------------------------------------------

# This script will trim the reads.

# Load modules  
fastp=/netfiles/thermofly/shared_software/fastp

#--------------------------------------------------------------------------------

WORKING_FOLDER=/gpfs2/scratch/jcnunez/thermofly/QC

#Define important file locations

meta=/users/j/c/jcnunez/scratch/thermofly/QC/Thermofly_metadata.tsv

#L=$(cat $meta | wc -l)
#echo ${L}

FIL1=$(cat ${meta} | awk -F '\t' '{print $24}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")
FIL2=$(cat ${meta} | awk -F '\t' '{print $25}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

SAMP_NAME=$(cat ${meta} | awk -F '\t' '{print $1}' |  sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d")

echo $FIL1
echo $FIL2
echo $SAMP_NAME

head -n1 $FIL1
head -n1 $FIL2

#--------------------------------------------------------------------------------
# This part of the pipeline will generate log files to record warnings and completion status

echo ${SAMP_NAME} "Trimming reads"
mkdir fastp_reports
mkdir cleaned_reads

# Call fastp and do some light trimming
$fastp \
-i ${FIL1} \
-I ${FIL2} \
-o $WORKING_FOLDER/cleaned_reads/${SAMP_NAME}_R1_clean.fq.gz \
-O $WORKING_FOLDER/cleaned_reads/${SAMP_NAME}_R2_clean.fq.gz \
--detect_adapter_for_pe \
--trim_front1 12 \
--trim_poly_g \
--thread 6 \
--cut_right \
--cut_right_window_size 6 \
--qualified_quality_phred 20 \
--html $WORKING_FOLDER/fastp_reports/${SAMP_NAME}_clean.html \
--json $WORKING_FOLDER/fastp_reports/${SAMP_NAME}_clean.json

# i = read 1
# I = read 2
# o & O = outputs
# outputs .json and html for qc

# Adapter trimming is enabled by default using overlap analysis, but for PE data you can also specify adapter sequence auto-detection by specifying --detect_adapter_for_pe
# For read 1 of PE data, the front trimming settings are --trim_front1
# Detect and trim polyG in read tails using --trim_poly_g
# Per read cutting by quality options:
### --cut_right: Move a sliding window from front to tail, if meet one window with mean quality < threshold, drop the bases in the window and the right part, then stop
### --cut_right_window_size: The window size option for cut_right
### --qualified_quality_phred: The quality vlue that a base is qualified 

#--------------------------------------------------------------------------------
# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo "pipeline completed" $(date)