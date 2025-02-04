# This script will annotate SNPs of interest

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

# Set working directory
setwd("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/Thermofly/6_Annotation")

# Load gene positions
SNPs <- read.csv("SNPs_of_interst.unq.txt", sep="\t")

# Load GDS file
gds.fn <- seqExampleFileName("D.basisetae.GATK.pipe.gds")
genofile.path <- "/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/Thermofly/6_Annotation/
genofile <- seqOpen(genofile.path)

###my GDS looks like:
# +    [  ] *
#   |--+ description   [  ] *
#   |--+ sample.id   { Str8 137 ZIP_ra(7.14%), 982B } *
#   |--+ variant.id   { Int32 10167176 ZIP_ra(34.6%), 13.4M } *
#   |--+ position   { Int32 10167176 ZIP_ra(40.7%), 15.8M } *
#   |--+ chromosome   { Str8 10167176 ZIP_ra(0.19%), 289.3K } *
#   |--+ allele   { Str8 10167176 ZIP_ra(16.4%), 6.4M } *
#   |--+ genotype   [  ] *
#   |  |--+ data   { Bit2 2x137x10167176 ZIP_ra(34.9%), 232.0M } *
#   |  |--+ extra.index   { Int32 3x0 ZIP_ra, 16B } *
#   |  \--+ extra   { Int16 0 ZIP_ra, 16B }
# |--+ phase   [  ]
# |  |--+ data   { Bit1 137x10167176 ZIP_ra(0.10%), 165.4K } *
#   |  |--+ extra.index   { Int32 3x0 ZIP_ra, 16B } *
#   |  \--+ extra   { Bit1 0 ZIP_ra, 16B }
# |--+ annotation   [  ]
# |  |--+ id   { Str8 10167176 ZIP_ra(0.10%), 9.7K } *
#   |  |--+ qual   { Float32 10167176 ZIP_ra(11.9%), 4.6M } *
#   |  |--+ filter   { Int32,factor 10167176 ZIP_ra(0.10%), 38.6K } *
#   |  |--+ info   [  ]
# |  |  |--+ INDEL   { Bit1 10167176 ZIP_ra(0.11%), 1.5K } *
#   |  |  |--+ IDV   { Int32 10167176 ZIP_ra(0.10%), 38.8K } *
#   |  |  |--+ IMF   { Float32 10167176 ZIP_ra(0.10%), 38.8K } *
#   |  |  |--+ DP   { Int32 10167176 ZIP_ra(41.5%), 16.1M } *
#   |  |  |--+ VDB   { Float32 10167176 ZIP_ra(93.5%), 36.3M } *
#   |  |  |--+ RPB   { Float32 10167176 ZIP_ra(91.2%), 35.4M } *
#   |  |  |--+ MQB   { Float32 10167176 ZIP_ra(73.7%), 28.6M } *
#   |  |  |--+ BQB   { Float32 10167176 ZIP_ra(86.7%), 33.6M } *
#   |  |  |--+ MQSB   { Float32 10167176 ZIP_ra(72.2%), 28.0M } *
#   |  |  |--+ SGB   { Float32 10167176 ZIP_ra(91.7%), 35.5M } *
#   |  |  |--+ MQ0F   { Float32 10167176 ZIP_ra(26.7%), 10.3M } *
#   |  |  |--+ ICB   { Float32 10167176 ZIP_ra(65.7%), 25.5M } *
#   |  |  |--+ HOB   { Float32 10167176 ZIP_ra(62.8%), 24.4M } *
#   |  |  |--+ AC   { Int32 10167176 ZIP_ra(34.7%), 13.4M } *
#   |  |  |--+ AN   { Int32 10167176 ZIP_ra(13.8%), 5.4M } *
#   |  |  |--+ DP4   { Int32 4x10167176 ZIP_ra(41.3%), 64.1M } *
#   |  |  |--+ MQ   { Int32 10167176 ZIP_ra(8.57%), 3.3M } *
#   |  |  |--+ ANN   { Str8 21554892 ZIP_ra(5.71%), 143.6M } *
#   |  |  |--+ LOF   { Str8 2667 ZIP_ra(19.6%), 16.2K } *
#   |  |  \--+ NMD   { Str8 1327 ZIP_ra(20.0%), 8.2K } *
#   |  \--+ format   [  ]
# |     \--+ PL   [  ] *
#   |        \--+ data   { VL_Int 137x30501528 ZIP_ra(55.9%), 2.9G } *
#   \--+ sample.annotation   [  ]

variant_ids <- seqGetData(genofile, "variant.id")
positions <- seqGetData(genofile, "position")
chromosomes <- seqGetData(genofile, "chromosome")

# Convert extracted data to a data frame
variant_info <- data.frame(variant_id = variant_ids, pos = positions, chr = chromosomes)

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



