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
data.binary.0.65.SNP <- read.table("Nucella_SNPs_maf0.05_pctind0.65_mindepth0.3_maxdepth2.binary.SNPs.gwas.lrt0", header = T, sep = "\t")
str(data.binary.0.65.SNP)

# Create unique Chromosome number 
Chr.unique <- unique(data.binary.0.65.SNP$Chromosome)
data.binary.0.65.SNP$CHR <- as.numeric(factor(data.binary.0.65.SNP$Chromosome, levels = Chr.unique))

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fails one of the filters) and are negative
data.binary.0.65.SNP.filt <- data.binary.0.65.SNP[-c(which(data.binary.0.65.SNP$LRT == -999), which(data.binary.0.65.SNP$LRT <= 0)), ]

hist(data.binary.0.65.SNP.filt$LRT, breaks = 50)

# ================================================================================== #

# Prepare data

# Name each SNP 
data.binary.0.65.SNP.filt$SNP <- paste("r", 1:length(data.binary.0.65.SNP.filt$Chromosome), sep="")

# Name each BP 
data.binary.0.65.SNP.filt$BP <- data.binary.0.65.SNP.filt$Position

# Get pvalues
data.binary.0.65.SNP.filt$P <- pchisq(data.binary.0.65.SNP.filt$LRT, df=1, lower=F)

# ================================================================================== #

# Make manhattan plot
manhattan(data.binary.0.65.SNP.filt, chr="CHR", bp="Position", p="P")

# Look at qq-plot of pvalues to check model fit
qqnorm(data.binary.0.65.SNP.filt$P)

# ================================================================================== #






