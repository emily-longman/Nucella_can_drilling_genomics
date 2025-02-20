# This script will annotate SNPs of interest.

# Prior to running this code, make sure the SnpEff vcf file from the previous step and the SNPs of interest are all in a directory called "Annotation"

#--------------------------------------------------------------------------------

# Load packages
install.packages(c('data.table', 'tidyverse'))
library(data.table)
library(tidyverse)
library(foreach)

# Load SeqArray
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

#--------------------------------------------------------------------------------

# Convert VCF to GDS

# Load the VCF file
vcf.fn <- "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/SNPeff/Nucella_SNPs.annotate.vcf"
# Parse the header
seqVCF_Header(vcf.fn)
# Convert VCF to GDS
seqVCF2GDS(vcf.fn, "Nucella_SNPs.annotate.gds")

#--------------------------------------------------------------------------------

# Open the GDS file
genofile <- seqOpen("/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/SNPeff/Nucella_SNPs.annotate.gds")
# Display the contents of the GDS file
genofile

# Structure of GDS file:
# +    [  ] *
#|--+ description   [  ] *
#|--+ sample.id   { Str8 192 LZMA_ra(5.27%), 1.2K } *
#|--+ variant.id   { Int32 168581 LZMA_ra(6.48%), 42.7K } *
#|--+ position   { Int32 168581 LZMA_ra(49.4%), 325.1K } *
#|--+ chromosome   { Str8 168581 LZMA_ra(1.09%), 23.2K } *
#|--+ allele   { Str8 168581 LZMA_ra(11.9%), 78.6K } *
#|--+ genotype   [  ] *
#|  |--+ data   { Bit2 2x192x168581 LZMA_ra(16.3%), 2.5M } *
#|  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
#|  \--+ extra   { Int16 0 LZMA_ra, 18B }
#|--+ phase   [  ]
#|  |--+ data   { Bit1 192x168581 LZMA_ra(0.02%), 737B } *
#|  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
#|  \--+ extra   { Bit1 0 LZMA_ra, 18B }
#|--+ annotation   [  ]
#|  |--+ id   { Str8 168581 LZMA_ra(0.10%), 181B } *
#|  |--+ qual   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |--+ filter   { Int32,factor 168581 LZMA_ra(0.04%), 253B } *
#|  |--+ info   [  ]
#|  |  |--+ INDEL   { Bit1 168581 LZMA_ra(0.65%), 145B } *
#|  |  |--+ IDV   { Int32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ IMF   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ DP   { Int32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ VDB   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ RPB   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ MQB   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ BQB   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ MQSB   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ NS   { Int32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ RPB2   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ MQB2   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ BQB2   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ MQSB2   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ SGB   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ MQ0F   { Float32 168581 LZMA_ra(0.04%), 253B } *
#|  |  |--+ I16   { Float32 16x168581 LZMA_ra(0.02%), 1.7K } *
#|  |  |--+ QS   { Float32 0 LZMA_ra, 18B } *
#|  |  |--+ AF   { Float32 168581 LZMA_ra(87.7%), 577.3K } *
#|  |  |--+ DPR   { Int32 0 LZMA_ra, 18B } *
#|  |  |--+ AD   { Int32 0 LZMA_ra, 18B } *
#|  |  |--+ ADF   { Int32 0 LZMA_ra, 18B } *
#|  |  |--+ ADR   { Int32 0 LZMA_ra, 18B } *
#|  |  |--+ ANN   { Str8 812106 LZMA_ra(6.50%), 7.9M } *
#|  |  |--+ LOF   { Str8 940 LZMA_ra(18.7%), 9.2K } *
#|  |  \--+ NMD   { Str8 509 LZMA_ra(20.1%), 5.3K } *
#|  \--+ format   [  ]
#|     |--+ DP4   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     |--+ SP   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     |--+ AD   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     |--+ ADF   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     |--+ ADR   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     |--+ PL   [  ] *
#|     |  \--+ data   { VL_Int 192x505743 LZMA_ra(18.3%), 19.3M } *
#|     |--+ DP   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     |--+ DV   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     |--+ DPR   [  ] *
#|     |  \--+ data   { VL_Int 192x0 LZMA_ra, 18B } *
#|     \--+ GL   [  ] *
#|        \--+ data   { Float32 192x505743 LZMA_ra(9.40%), 34.8M } *
#\--+ sample.annotation   [  ]

#--------------------------------------------------------------------------------

# Load SNPs of interest
GWAS_snps <- read.csv("Nucella_GWAS_outlier_SNPs.csv", header=T)

# Create SNP ID column
snps <- GWAS_snps %>% mutate(SNP_id = paste(Chromosome, Position, sep = "_"))

str(snps)

#--------------------------------------------------------------------------------

# Extract SNP data from GDS
snp.dt <- data.table(Chromosome=seqGetData(genofile, "chromosome"),
                        Position=seqGetData(genofile, "position"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        id=seqGetData(genofile, "variant.id")) %>%
                        mutate(SNP_id = paste(Chromosome, Position, sep = "_"))


#--------------------------------------------------------------------------------


# Extract annotation data for each SNP of interest

annotation = 
foreach(i=1:dim(snps)[1], 
.combine = "rbind",
.errorhandling = "remove")%do%{

message(i)
seqResetFilter(genofile)

tmp.i = snps[i,]$SNP_id

pos.tmp = snp.dt %>% filter(SNP_id %in% tmp.i) %>% .$id

#if (identical(pos.tmp, integer(0))) {print("NA")

#} else {

seqSetFilter(genofile, variant.id = pos.tmp)

ann_data <- seqGetData(genofile, "annotation/info/ANN")$data

L = length(ann_data)

annotate.list =
foreach(k=1:L, 
.combine = "rbind")%do%{

  tmp = ann_data[k] 
  tmp2= str_split(tmp, "\\|")
  
  data.frame(
    id=pos.tmp,
    SNP_id = tmp.i,
    annotation.id=k,
    Allele = tmp2[[1]][1],
    Annotation = tmp2[[1]][2],
    Annotation_Impact = tmp2[[1]][3],
    Gene_Name = tmp2[[1]][4],
    Gene_ID = tmp2[[1]][5],
    Feature_Type = tmp2[[1]][6],
    Feature_ID = tmp2[[1]][7],
    Transcript_BioType = tmp2[[1]][8],
    Rank = tmp2[[1]][9],
    HGVS.c = tmp2[[1]][10],
    HGVS.p = tmp2[[1]][11],
    cDNA.pos.cDNA.length = tmp2[[1]][12],
    CDS.pos.CDS.length = tmp2[[1]][13],
    AA.pos.AA.length = tmp2[[1]][14],
    Distance = tmp2[[1]][15]

  )
}

return(annotate.list)
#}
}



#--------------------------------------------------------------------------------

# Join annotation and SNP information
N.canaliculata_annotated_SNPs <- left_join(annotation, snp.dt, by = join_by(SNP_id))

# Write output
write.csv(N.canaliculata_annotated_SNPs, "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/SNPeff/N.canaliculata_SNPs_annotated.txt")
write.csv(N.canaliculata_annotated_SNPs, "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/fastq_to_GL/SNPeff/N.canaliculata_SNPs_annotated.csv")









#####

#example from D. bas to search genome for specific nucleotide sequences.
samtools faidx GCA_035041595.1_ASM3504159v1_genomic.fna.masked.fa "JAWNLB010000001.1":5-6