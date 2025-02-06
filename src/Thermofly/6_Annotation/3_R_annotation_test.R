# This script will annotate SNPs of interest

# Create Annotation dir in working_directory to test things out
working_folder=/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly
cd $working_folder
mkdir annotation
scp SnpEff/D.bas_SNPs.annotate.vcf annotation/
cd annotation

#--------------------------------------------------------------------------------

module load R/4.4.1  

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


#--------------------------------------------------------------------------------

# Convert VCF to GDS

# Load the VCF file
vcf.fn <- "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly/annotation/D.bas_SNPs.annotate.vcf"
# Parse the header
seqVCF_Header(vcf.fn)
# Convert VCF to GDS
seqVCF2GDS(vcf.fn, "D.bas_SNPs.annotate.gds")

#--------------------------------------------------------------------------------

# Load GDS file
#gds.fn <- seqExampleFileName("gds")
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


#--------------------------------------------------------------------------------

# Extract annotation from GDS
#ann_data <- seqGetData(geno.file, "annotation/info/ANN")



#--------------------------------------------------------------------------------


# Load SNPs of interest
snps <- read.csv("SNPs_of_interst.unq.txt", sep="\t")
str(snps)

colnames(snps)[colnames(snps) == 'pos.id'] <- 'id'
#snps$chr.pos <-  paste(snps$chr, ".", snps$pos)

#--------------------------------------------------------------------------------

# Extract SNP data from GDS
snp.dt <- data.table(chr=seqGetData(geno.file, "chromosome"),
                        pos=seqGetData(geno.file, "position"),
                        nAlleles=seqGetData(geno.file, "$num_allele"),
                        id=seqGetData(geno.file, "variant.id"))

# Create id that combines chromosome and position
#snp.dt$SNP_id <- paste0(snp.dt$chr, "_", snp.dt$pos)

#--------------------------------------------------------------------------------

# Extract annotation data for each SNP of interest

annotation = 
foreach(i=1:dim(snps)[1], .combine = "rbind")%do%{
seqResetFilter(geno.file)

myvar = snp.dt %>% filter(chr == CHR & pos == POS)$id

seqSetFilter(geno.file, variant.id= myvar )
ann_data <- seqGetData(geno.file, "annotation/info/ANN")$data

L = length(ann_data)
foreach(k=1:L, .combine = "rbind")%do%{
  tmp = ann_data[k] 
  tmp2= str_split(tmp, "\\|")
  

  data.frame(
    id=i,
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


ERRORS_WARNINGS_INFO = tmp2[[1]][16],

# Join annotation and SNP information
D.basisetae_annotated_SNPs <- left_join(annotation, snp.dt, by = join_by(id))


write.csv(D.basisetae_annotated_SNPs, "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly/annotation/D.basisetae_annotated_SNPs.txt")

write.csv(D.basisetae_annotated_SNPs, "/gpfs2/scratch/elongman/Nucella_can_drilling_genomics/data/processed/Thermofly/annotation/D.basisetae_annotated_SNPs.csv")

#--------------------------------------------------------------------------------






# Example code:

# SNPs of interest
adapt_cands = snp.dt[sample(10),]

ann_ex = 
foreach(i=1:dim(adapt_cands)[1], .combine = "rbind")%do%{
seqResetFilter(geno.file)

myvar = adapt_cands$id[i]

seqSetFilter(geno.file, variant.id= myvar )
ann_data <- seqGetData(geno.file, "annotation/info/ANN")$data

L = length(ann_data)
foreach(k=1:L, .combine = "rbind")%do%{
  tmp = ann_data[k] 
  tmp2= str_split(tmp, "\\|")
  

  data.frame(
    id=i,
    annotation.id=k,
    allele=tmp2[[1]][1],
    feature = tmp2[[1]][2]

  )

}
}

left_join(ann_ex, snp.dt)





#--------------------------------------------------------------------------------
# Extract data from GDS
variant_ids <- seqGetData(geno.file, "variant.id")
positions <- seqGetData(geno.file, "position")
chromosomes <- seqGetData(geno.file, "chromosome")

# Convert extracted data to a data frame
variant_info <- data.frame(variant_id = variant_ids, pos = positions, chr = chromosomes)



# Separate the 'data' and 'length' components
ann_entries <- ann_data$data
ann_lengths <- ann_data$length

# Initialize a list placeholder to store the extracted data for each variant
ann_list <- vector("list", length = length(ann_lengths))

# Set an index to keep track of position in `ann_entries`
index <- 1

# Loop through each variant using the length of each entry to slice ann_entries
for (i in seq_along(ann_lengths)) {
  message("Processing variant ID: ", i)
  
  # Number of entries for this variant
  n_entries <- ann_lengths[i]
  
  # Extract the entries for this variant
  if (n_entries > 0) {
    variant_ann <- ann_entries[index:(index + n_entries - 1)]
    ann_list[[i]] <- strsplit(variant_ann, "\\|")  # Split each annotation by '|'
    index <- index + n_entries
  } else {
    ann_list[[i]] <- list(rep(NA, 15))  # If no annotation, assign NA
  }
}

# Flatten and convert the list into a data frame
ann_df <- do.call(rbind, lapply(ann_list, function(x) {
  do.call(rbind, lapply(x, function(y) if (!is.na(y)) y else rep(NA, 15)))
}))


# Set column names for the ann_df data frame
colnames(ann_df) <- c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", 
                      "Gene_ID", "Feature_Type", "Feature_ID", 
                      "Transcript_BioType", "Rank", "HGVS.c", 
                      "HGVS.p", "cDNA.pos", "CDS.pos", 
                      "AA.pos", "Distance", "ERRORS_WARNINGS_INFO")

# Convert ann_df to a data frame
ann_df <- as.data.frame(ann_df, stringsAsFactors = FALSE)

# Ensure variant_info and ann_df have the same number of rows before combining
if (nrow(variant_info) == nrow(ann_df)) {
  # Combine variant information with annotation data
  final_df <- cbind(variant_info, ann_df)
} else {
  stop("Number of variants in variant_info does not match number of rows in ann_df")
}

# Close the GDS file
seqClose(genofile)

#--------------------------------------------------------------------------------

# Check the combined data
head(final_df)

###load SNP list interested in
significant_windows<-load("windowed.RData")


###mine looks like this:
#test_code            chr           variable pos_mean pos_min pos_max win
#1         0 NW_022145594.1 ALLSKY_SFC_UVB_min    54798    7823  102218   1
#2         1 NW_022145594.1 ALLSKY_SFC_UVB_min    54798    7823  102218   1
#3         2 NW_022145594.1 ALLSKY_SFC_UVB_min    54798    7823  102218   1
#rnp.pr.0.05 rnp.pr.0.01 rnp.pr.0.001 rnp.binom.p.0.05 rnp.binom.p.0.01
#1  0.01620370 0.001157407  0.000000000     2.881103e-07      0.003111481
#2  0.03587963 0.009259259  0.001157407     6.045702e-02      1.000000000
#3  0.04166667 0.005787037  0.000000000     3.093977e-01      0.300537578

###used fuzzy join to join my windowed analysis (simpler if its 1:1 if looking at individual SNPS)
filtered_snps <- final_df %>%
  fuzzy_inner_join(
    significant_windows,
    by = c("chr" = "chr", "variable" = "variable", "pos" = "pos_min", "pos" = "pos_max"),
    match_fun = list(`==`, `==`, `>=`, `<=`)
  ) %>%
  distinct()
