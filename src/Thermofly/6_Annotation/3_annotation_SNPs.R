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
vcf.fn <- "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly/annotation/D.bas_SNPs.annotate.vcf"
# Parse the header
seqVCF_Header(vcf.fn)
# Convert VCF to GDS
seqVCF2GDS(vcf.fn, "D.bas_SNPs.annotate.gds")

#--------------------------------------------------------------------------------

# Open the GDS file
geno.file <- seqOpen("/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly/annotation/D.bas_SNPs.annotate.gds")
# Display the contents of the GDS file
geno.file

# Structure of GDS file:
# +    [  ] *
#|--+ description   [  ] *
#|--+ sample.id   { Str8 18 LZMA_ra(22.9%), 189B } *
#|--+ variant.id   { Int32 4520251 LZMA_ra(3.00%), 530.4K } *
#|--+ position   { Int32 4520251 LZMA_ra(24.9%), 4.3M } *
#|--+ chromosome   { Str8 4520251 LZMA_ra(0.02%), 12.2K } *
#|--+ allele   { Str8 4520251 LZMA_ra(23.1%), 6.0M } *
#|--+ genotype   [  ] *
#|  |--+ data   { Bit2 2x18x4582620 LZMA_ra(30.6%), 12.0M } *
#|  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
#|  \--+ extra   { Int16 0 LZMA_ra, 18B }
#|--+ phase   [  ]
#|  |--+ data   { Bit1 18x4520251 LZMA_ra(32.8%), 3.2M } *
#|  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
#|  \--+ extra   { Bit1 0 LZMA_ra, 18B }
#|--+ annotation   [  ]
#|  |--+ id   { Str8 4520251 LZMA_ra(0.02%), 813B } *
#|  |--+ qual   { Float32 4520251 LZMA_ra(65.3%), 11.3M } *
#|  |--+ filter   { Int32,factor 4520251 LZMA_ra(0.02%), 2.7K } *
#|  |--+ info   [  ]
#|  |  |--+ AC   { Int32 5046291 LZMA_ra(10.9%), 2.1M } *
#|  |  |--+ AF   { Float32 5046291 LZMA_ra(18.7%), 3.6M } *
#|  |  |--+ AN   { Int32 4520251 LZMA_ra(8.05%), 1.4M } *
#|  |  |--+ BaseQRankSum   { Float32 4520251 LZMA_ra(32.4%), 5.6M } *
#|  |  |--+ DP   { Int32 4520251 LZMA_ra(22.1%), 3.8M } *
#|  |  |--+ END   { Int32 4520251 LZMA_ra(0.02%), 2.7K } *
#|  |  |--+ ExcessHet   { Float32 4520251 LZMA_ra(14.2%), 2.5M } *
#|  |  |--+ FS   { Float32 4520251 LZMA_ra(26.6%), 4.6M } *
#|  |  |--+ InbreedingCoeff   { Float32 4520251 LZMA_ra(46.8%), 8.1M } *
#|  |  |--+ MLEAC   { Int32 5046291 LZMA_ra(11.0%), 2.1M } *
#|  |  |--+ MLEAF   { Float32 5046291 LZMA_ra(18.5%), 3.6M } *
#|  |  |--+ MQ   { Float32 4520251 LZMA_ra(4.07%), 717.8K } *
#|  |  |--+ MQRankSum   { Float32 4520251 LZMA_ra(4.25%), 751.0K } *
#|  |  |--+ QD   { Float32 4520251 LZMA_ra(41.4%), 7.1M } *
#|  |  |--+ RAW_MQandDP   { Int32 2x4520251 LZMA_ra(0.01%), 5.3K } *
#|  |  |--+ ReadPosRankSum   { Float32 4520251 LZMA_ra(35.0%), 6.0M } *
#|  |  |--+ SOR   { Float32 4520251 LZMA_ra(33.6%), 5.8M } *
#|  |  |--+ ANN   { Str8 7957439 LZMA_ra(6.50%), 52.7M } *
#|  |  |--+ LOF   { Str8 10017 LZMA_ra(3.28%), 6.7K } *
#|  |  \--+ NMD   { Str8 1425 LZMA_ra(7.36%), 2.1K } *
#|  \--+ format   [  ]
#|     |--+ AD   [  ] *
#|     |  \--+ data   { VL_Int 18x9566542 LZMA_ra(31.6%), 54.6M } *
#|     |--+ DP   [  ] *
#|     |  \--+ data   { VL_Int 18x4520251 LZMA_ra(40.9%), 33.1M } *
#|     |--+ GQ   [  ] *
#|     |  \--+ data   { VL_Int 18x4520251 LZMA_ra(47.2%), 42.0M } *
#|     |--+ MIN_DP   [  ] *
#|     |  \--+ data   { VL_Int 18x0 LZMA_ra, 18B } *
#|     |--+ PGT   [  ] *
#|     |  \--+ data   { Str8 18x2282075 LZMA_ra(6.85%), 4.0M } *
#|     |--+ PID   [  ] *
#|     |  \--+ data   { Str8 18x2282075 LZMA_ra(8.40%), 10.0M } *
#|     |--+ PL   [  ] *
#|     |  \--+ data   { VL_Int 18x15295094 LZMA_ra(45.3%), 169.1M } *
#|     |--+ PS   [  ] *
#|     |  \--+ data   { VL_Int 18x2282075 LZMA_ra(4.72%), 8.9M } *
#|     |--+ RGQ   [  ] *
#|     |  \--+ data   { VL_Int 18x0 LZMA_ra, 18B } *
#|     \--+ SB   [  ] *
#|        \--+ data   { VL_Int 18x0 LZMA_ra, 18B } *
#\--+ sample.annotation   [  ]

#--------------------------------------------------------------------------------

# Load SNPs of interest
snps <- read.csv("FSTGLM_SNPS.outlier.Feb7.V2.txt", sep="\t")
str(snps)

#--------------------------------------------------------------------------------

# Extract SNP data from GDS
snp.dt <- data.table(chr=seqGetData(geno.file, "chromosome"),
                        pos=seqGetData(geno.file, "position"),
                        nAlleles=seqGetData(geno.file, "$num_allele"),
                        id=seqGetData(geno.file, "variant.id")) %>%
                        mutate(SNP_id = paste(chr, pos, sep = "_"))


#--------------------------------------------------------------------------------

# Extract annotation data for each SNP of interest

annotation = 
foreach(i=1:dim(snps)[1], 
.combine = "rbind")%do%{

seqResetFilter(geno.file)

tmp.i = snps[i,]$SNP_id

pos.tmp = snp.dt %>% filter(SNP_id %in% tmp.i) %>% .$id

seqSetFilter(geno.file, variant.id = pos.tmp)

ann_data <- seqGetData(geno.file, "annotation/info/ANN")$data

L = length(ann_data)
foreach(k=1:L, .combine = "rbind")%do%{
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
}

#--------------------------------------------------------------------------------

# Join annotation and SNP information (can join by "id" or "SNP_ic")
D.basisetae_annotated_SNPs <- left_join(annotation, snp.dt, by = join_by(id))


# Write output
write.csv(D.basisetae_annotated_SNPs, "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly/annotation/D.basisetae_annotated_SNPs.txt")
write.csv(D.basisetae_annotated_SNPs, "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly/annotation/D.basisetae_annotated_SNPs.csv")