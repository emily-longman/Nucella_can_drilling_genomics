## Plotting with plot_admix_function.R code

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Note: to run script move gwas output to results dir 

# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

# Set relative path of results directory from root
dir(find_root_file("results", criterion = has_file("README.md")))
results_path_from_root <- find_root_file("results", "stats", "gwas", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(results_path_from_root)

# Set working directory as path from root
setwd(results_path_from_root)

# ================================================================================== #

# This script will make a manhattan plot. 

# ================================================================================== #

# Load packages
library(qqman)

# ================================================================================== #

# Load data 
data.binary.SNP.minCount <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.binary.SNPs.minCount20.gwas.lrt0", header = T, sep = "\t")
str(data.binary.SNP.minCount)

# Create unique Chromosome number 
Chr.unique <- unique(data.binary.SNP.minCount$Chromosome)
data.binary.SNP.minCount$CHR <- as.numeric(factor(data.binary.SNP.minCount$Chromosome, levels = Chr.unique))

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fails one of the filters) and are negative
data.binary.SNP.minCount.filt <- data.binary.SNP.minCount[-c(which(data.binary.SNP.minCount$LRT == -999), which(data.binary.SNP.minCount$LRT <= 0)), ]

hist(data.binary.SNP.minCount.filt$LRT, breaks = 50)

# ================================================================================== #

# Prepare data

# Name each SNP 
data.binary.SNP.minCount.filt$SNP <- paste("r", 1:length(data.binary.SNP.minCount.filt$Chromosome), sep="")

# Name each BP 
data.binary.SNP.minCount.filt$BP <- data.binary.SNP.minCount.filt$Position

# Get pvalues
data.binary.SNP.minCount.filt$P <- pchisq(data.binary.SNP.minCount.filt$LRT, df=1, lower=F)

# ================================================================================== #

# Make manhattan plot
manhattan(data.binary.SNP.minCount.filt, chr="CHR", bp="Position", p="P")

# Look at qq-plot of pvalues to check model fit
qqnorm(data.binary.SNP.minCount.filt$P)

# ================================================================================== #

# Make a list of the candidate loci
candidates <- data.binary.SNP.minCount.filt[which(data.binary.SNP.minCount.filt$LRT > 15),]$SNP


manhattan(data.binary.SNP.minCount.filt, chr="CHR", bp="Position", p="P", highlight=candidates,  cex=0.6)


