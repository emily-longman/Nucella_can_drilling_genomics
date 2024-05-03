#J loop#	#I will now extract some summary stats
samtools flagstat --threads $CPU \ 
$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
> $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt

#J loop#	#build bam files
samtools view -b -q $QUAL --threads $CPU  \ 
$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam \
> $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam

#J loop#	# Sort with picard
# Notice that once a file has been sorted it is added the "srt" suffix
java -Xmx$JAVAMEM -jar $PICARD SortSam \ 
I=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam \
O=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam \
SO=coordinate \
VALIDATION_STRINGENCY=SILENT

#J loop# Remove duplicates with picard
	# Notice that once a file has been sorted it is added the "rmdp" suffix
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \ 
I=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam \
O=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.rmdp.bam \
M=$WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.dupstat.txt \
VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#J loop# Lets do QC on the bam file
qualimap bamqc -bam $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.rmdp.bam -outdir $WORKING_FOLDER/mapping_stats/Qualimap_${i} \ 
--java-mem-size=$JAVAMEM

#J loop#	# Clean intermediate files
rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.sam
rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.bam
rm $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.srt.bam

#J loop#	# Housekeeping
mv $WORKING_FOLDER/${j}_reads/${i}/${i}.flagstats_raw_${j}.sam.txt \ 
$WORKING_FOLDER/mapping_stats
mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.dupstat.txt \ 
$WORKING_FOLDER/mapping_stats

mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim_fastqc.html \ 
$WORKING_FOLDER/read_stats
mv $WORKING_FOLDER/${j}_reads/${i}/${i}.${j}.reads.strict.trim_fastqc.zip \ 
$WORKING_FOLDER/read_stats

#J loop#	
done # End loop of j

#--------------------------------------------------------------------------------

# Merge and assess the final file

# Here I will merge the bam outputs from the merge and unmerged portions of the pipeline. 
# Subsequently, I will once again sort and remove duplicated, before performing the final QC on the aligment.

# Merge bams
java -Xmx$JAVAMEM -jar $PICARD MergeSamFiles \ 
I=$WORKING_FOLDER/merged_reads/${i}/${i}.merged.srt.rmdp.bam  \ 
I=$WORKING_FOLDER/unmerged_reads/${i}/${i}.unmerged.srt.rmdp.bam  \ 
O=$WORKING_FOLDER/joint_bams/${i}.joint.bam

# Sort merge bams
java -Xmx$JAVAMEM -jar $PICARD SortSam \ 
I=$WORKING_FOLDER/joint_bams/${i}.joint.bam \ 
O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \ 
SO=coordinate \ 
VALIDATION_STRINGENCY=SILENT

# Remove duplicates of final file
java -Xmx$JAVAMEM -jar $PICARD MarkDuplicates \ 
I=$WORKING_FOLDER/joint_bams/${i}.joint.srt.bam \ 
O=$WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam  \ 
M=$WORKING_FOLDER/mapping_stats/${i}.joint.dupstat.txt \ 
VALIDATION_STRINGENCY=SILENT \ 
REMOVE_DUPLICATES=true

# Assess quality of final file
qualimap bamqc -bam $WORKING_FOLDER/joint_bams/${i}.joint.srt.rmdp.bam  -outdir $WORKING_FOLDER/joint_bams_qualimap/Qualimap_JointBam_${i} \ 
--java-mem-size=$JAVAMEM
 
# Remove intermediary files
rm $WORKING_FOLDER/joint_bams/${i}.joint.bam
rm $WORKING_FOLDER/joint_bams/${i}.joint.srt.bam

###########################################################################
###########################################################################
