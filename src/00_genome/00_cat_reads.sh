#!/bin/bash

# Concatenate reads; all pass filter

# Move to the directory where you want the output files to be saved
cd /netfiles/pespenilab_share/Nucella/raw/ONT/

# Concatenate reads from each ONT cell
cat /netfiles/pespenilab_share/Nucella/raw/ONT/PROM0150_Sanford_NC3-PCR_FC2_11152023/20231115_1701_2E_PAS52172_fe64f551/fastq_pass/*.fastq.gz > FC1.ONT.nuc.fastq.gz

cat /netfiles/pespenilab_share/Nucella/raw/ONT/PROM0150_Sanford_NC3-PCR_05312023/20230531_1740_3G_PAQ05812_93506d2a/fastq_pass/*.fastq.gz > FC2.ONT.nuc.fastq.gz

# Concatenate reads into one fastq.gz file
cat FC1.ONT.nuc.fastq.gz FC2.ONT.nuc.fastq.gz > FC_all.ONT.nuc.fastq.gz
