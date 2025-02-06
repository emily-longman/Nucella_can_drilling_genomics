# This script will annotate SNPs of interest

#--------------------------------------------------------------------------------

# Load packages
install.packages(c('data.table', 'tidyverse'))
library(data.table)
library(tidyverse)
#library(GWASTools)

# Load SeqArray
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

#--------------------------------------------------------------------------------

# Set working directory
setwd("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/Thermofly/6_Annotation")

#--------------------------------------------------------------------------------

# Load gene positions
SNPs <- read.csv("SNPs_of_interst.unq.txt", sep="\t")

#--------------------------------------------------------------------------------

# Convert VCF to GDS

# Load the VCF file
vcf.fn <- "/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/Thermofly/6_Annotation/D.bas_SNPs.annotate.vcf"
# Parse the header
seqVCF_Header(vcf.fn)
# Convert VCF to GDS
seqVCF2GDS(vcf.fn, "D.bas_SNPs.annotate.gds")

#--------------------------------------------------------------------------------

# Load GDS file
gds.fn <- seqExampleFileName("gds")
# Open the GDS file
geno.file <- seqOpen(gds.fn)
# Display the contents of the GDS file
geno.file

# Structure of GDS file:
# +    [  ] *
#|--+ description   [  ] *
#  |--+ sample.id   { Str8 90 LZMA_ra(34.7%), 257B } *
#  |--+ variant.id   { Int32 1348 LZMA_ra(16.7%), 905B } *
#  |--+ position   { Int32 1348 LZMA_ra(64.4%), 3.4K } *
#  |--+ chromosome   { Str8 1348 LZMA_ra(4.39%), 157B } *
#  |--+ allele   { Str8 1348 LZMA_ra(16.6%), 901B } *
#  |--+ genotype   [  ] *
#  |  |--+ data   { Bit2 2x90x1348 LZMA_ra(26.3%), 15.6K } *
#  |  |--+ ~data   { Bit2 2x1348x90 LZMA_ra(29.2%), 17.3K } *
#  |  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
#  |  \--+ extra   { Int16 0 LZMA_ra, 18B }
#|--+ phase   [  ]
#|  |--+ data   { Bit1 90x1348 LZMA_ra(0.86%), 137B } *
#  |  |--+ ~data   { Bit1 1348x90 LZMA_ra(0.86%), 137B } *
#  |  |--+ extra.index   { Int32 3x0 LZMA_ra, 18B } *
#  |  \--+ extra   { Bit1 0 LZMA_ra, 18B }
#|--+ annotation   [  ]
#|  |--+ id   { Str8 1348 LZMA_ra(38.3%), 5.5K } *
#  |  |--+ qual   { Float32 1348 LZMA_ra(2.11%), 121B } *
#  |  |--+ filter   { Int32,factor 1348 LZMA_ra(2.11%), 121B } *
#  |  |--+ info   [  ]
#|  |  |--+ AA   { Str8 1328 LZMA_ra(22.1%), 593B } *
#  |  |  |--+ AC   { Int32 1348 LZMA_ra(24.1%), 1.3K } *
#  |  |  |--+ AN   { Int32 1348 LZMA_ra(19.6%), 1.0K } *
#  |  |  |--+ DP   { Int32 1348 LZMA_ra(47.7%), 2.5K } *
#  |  |  |--+ HM2   { Bit1 1348 LZMA_ra(145.6%), 253B } *
#  |  |  |--+ HM3   { Bit1 1348 LZMA_ra(145.6%), 253B } *
#  |  |  |--+ OR   { Str8 1348 LZMA_ra(19.6%), 341B } *
#  |  |  |--+ GP   { Str8 1348 LZMA_ra(24.3%), 3.8K } *
#  |  |  \--+ BN   { Int32 1348 LZMA_ra(20.7%), 1.1K } *
#  |  \--+ format   [  ]
#|     \--+ DP   [  ] *
#  |        |--+ data   { VL_Int 90x1348 LZMA_ra(70.8%), 115.2K } *
#  |        \--+ ~data   { VL_Int 1348x90 LZMA_ra(65.1%), 105.9K } *
#  \--+ sample.annotation   [  ]
#\--+ family   { Str8 90 LZMA_ra(55.0%), 221B } *


#--------------------------------------------------------------------------------

# Extract data from GDS
variant_ids <- seqGetData(geno.file, "variant.id")
positions <- seqGetData(geno.file, "position")
chromosomes <- seqGetData(geno.file, "chromosome")
annot <- 
  
  seqGetData(geno.file, "annotation/DP/data")

message("Annotations")
tmp <- seqGetData(geno.file, "annotation/info/ANN")



# Close the GDS file
#seqClose(geno.file)

#--------------------------------------------------------------------------------


VCF.R
vcf@info$ANN


# Flatten and convert ann_list into a data frame
ann_df <- do.call(rbind, lapply(ann_list, function(x) {
  do.call(rbind, lapply(x, function(y) {
    if (length(y) > 15) {
      y[1:15]
    } else if (length(y) < 15) {
      c(y, rep(NA, 15 - length(y)))
    } else {
      y
    }
  }))
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



