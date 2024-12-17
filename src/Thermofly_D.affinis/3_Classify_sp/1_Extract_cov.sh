# Loop for mining coverage data from BAMQC files

# Loop to extract raw data

# Extract coverage data for each contig
cat genome_results.txt | sed -n '129,3862p;3862q' | sed 's/^\t//g' > test1.txt

# Add column with sample name
awk 'BEGIN{OFS="\t"} {print $0, (FNR>1 ? FILENAME : FILENAME)}' test1.txt | head

######

#Loop through Bam qc to create masker 

testfile=/netfiles/thermofly/ANALYSES/Affinis_Altitude/bamQCs/bams_qualimap/Qualimap_D_aff.wild.US-VT-Cam.1_7_2022.w.CH800D.1/genome_results.txt

cat $testfile | sed -n '129,3862p;3862q' | sed -e 's/^\t//g' > test1.txt

# Added column for each row that says "test1.txt"
awk 'BEGIN{OFS="\t"} {print $0, (FNR>1 ? FILENAME : FILENAME)}' test1.txt | head


# Loop through qualimap files
files=echo $(ls /netfiles/thermofly/ANALYSES/Affinis_Altitude/bamQCs/bams_qualimap/Qualimap_D_aff.wild.*/genome_results.txt)

files=/netfiles/thermofly/ANALYSES/Affinis_Altitude/bamQCs/bams_qualimap/Qualimap_D_aff.wild.*/genome_results.txt
for i in files
do
echo $i
cat $i | sed -n '129,3862p;3862q' | sed -e 's/^\t//g' > $i.out.text
awk 'BEGIN{OFS="\t"} {print $0, (FNR>1 ? FILENAME : FILENAME)}' $i.out.txt >> bam_extract.txt
done

# Open file and use find/replace to get rid of path so just name of sample


# Open in R
library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(FactoMineR)

# Load data
in_data <- fread("bam_extract.txt")
names(in_dat) = c("chr", "L", "mappedReads", "mean_cov", "sd_cov", "sample")

aff=c("JAZHFB)10000004.1", ..... etc.)

in_data %<>% mutate(sp_chrs = case_when(chr %in% aff ~ "Affinis",))

grep("JAWNKY", in_data$chr) -> alg_chr


# large chrom exploration
in_dat %>% filter(L > 100000) # most of Athabasca genome is fragmented so don't only the biggest

# cast long for PCA
in_dat %>% dcast(sample~chr, value.var = "mean_cov") -> in_dat_cast

in_dat_cast[-1] -> in_dat_cast.redux

# Add samples to rownames
rownames(in_dat_cast) = in_dat_cast[-1]

in_dat_cast.redux %>% PCA(graph = F) -> PCA.calc.obj

PCA.calc.obj$ind$coord %>% as.data.frame %>% mutate(sampleID = rownames(.)) -> PCA.coords

PCA.coords %>%
ggplot(aes(x=Dim.1, y=Dim.2)) + geom_point(size=3) -> firstpass.PCA.plot


# Manually check the chromosome names for each genome (for example for algonquin)
grep "^>" GCA_035041765.1_ASM3504176v1_genomic.fna.masked


