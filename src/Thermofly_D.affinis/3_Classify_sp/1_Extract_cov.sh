# Loop for mining coverage data from BAMQC files

# Loop to extract raw data

129-3862

cat genome_results.txt | sed -n '129,3862p;3862q' | sed 's/^\t//g' > test.1

awk -i 'BEGIN{OFS="\t} {print $0} ()