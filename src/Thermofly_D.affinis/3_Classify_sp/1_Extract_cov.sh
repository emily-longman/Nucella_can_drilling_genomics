# Loop for mining coverage data from BAMQC files

# Loop to extract raw data

129-3862

# Extract coverage data for each contig
cat genome_results.txt | sed -n '129,3862p;3862q' | sed 's/^\t//g' > test1.txt

# Add column with sample name


######

#JCBN

testfile=/netfiles/thermofly/ANALYSES/Affinis_Altitude/bamQCs/bams_qualimap/Qualimap_D_aff.wild.US-VT-Cam.1_7_2022.w.CH800D.1/genome_results.txt

cat $testfile | sed -n '129,3862p;3862q' | sed -e 's/^\t//g' > test1.txt

# Added column for each row that says "test1.txt"
awk 'BEGIN{OFS="\t"} {print $0, (FNR>1 ? FILENAME : FILENAME)}' test1.txt | head