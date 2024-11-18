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

# Get R manhattan plot function (https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function)
source("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/01_Fastq_to_GL/09_GWAS/02_gwas_graphing.R")

# ================================================================================== #

# Load data 
data.binary <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.binary.gwas.lrt0", header = T, sep = "\t")
str(data.binary)

# Create unique Chromosome number 
Chr.unique <- unique(data.binary$Chromosome)
data.binary$CHR <- as.numeric(factor(data.binary$Chromosome, levels = Chr.unique))

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fails one of the filters) and are negative
data.binary.filt <- data.binary[-c(which(data.binary$LRT == -999), which(data.binary$LRT <= 0)), ]

hist(data.binary.filt$LRT, breaks = 50)

# ================================================================================== #

# Prepare data

# Name each SNP 
data.binary.filt$SNP <- paste("r", 1:length(data.binary.filt$Chromosome), sep="")

# Name each BP 
data.binary.filt$BP <- data.binary.filt$Position

# Get pvalues
data.binary.filt$P <- pchisq(data.binary.filt$LRT, df=1, lower=F)

# Filter data (none of these in my data)
#data.binary.filt <- data.binary.filt[-c(which(data.binary.filt$pvalue == "NaN" ),
#                                        which(data.binary.filt$pvalue == "Inf"),
#                                        which(data.binary.filt$LRT == "Inf")),]

# ================================================================================== #

# Make manhattan plot
manhattan(data.binary.filt, chr="CHR", bp="Position", p="P", highlight = hlight)

# Look at qq-plot of pvalues to check model fit
qqnorm(data.binary.filt$P)

# ================================================================================== #

# Make manhattan plot prettier (https://github.com/pcgoddard/Burchardlab_Tutorials/wiki/GGplot2-Manhattan-Plot-Function)
# Not working - maybe only ok for chromosomes 

mypalette <- c("#5D82BB", "#3B64A5", "#1E4F9E", "#103B7E", "#082B64")
hlight <- data.binary.filt$SNP[which(data.binary.filt$P < 0.00003)]
sig = 5e-5 # significant threshold line
sugg = 1e-6 

#gg.manhattan(data.binary.filt, threshold=1e-6, hlight=hlight, col=NA, ylims=c(0,10), title="My Manhattan Plot")







