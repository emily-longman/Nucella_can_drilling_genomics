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
library(qqman) #https://cran.r-project.org/web/packages/qqman/vignettes/qqman.html

# ================================================================================== #

# Load data 
data.count.SNP <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.count.SNPs.gwas.lrt0", header = T, sep = "\t")
str(data.count.SNP)

# Create unique Chromosome number 
Chr.unique <- unique(data.count.SNP$Chromosome)
data.count.SNP$CHR <- as.numeric(factor(data.count.SNP$Chromosome, levels = Chr.unique))

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fails one of the filters) and are negative
data.count.SNP.filt <- data.count.SNP[-c(which(data.count.SNP$LRT == -999), which(data.count.SNP$LRT <= 0)), ]

hist(data.count.SNP.filt$LRT, breaks = 50)

# ================================================================================== #

# Prepare data

# Name each SNP 
data.count.SNP.filt$SNP <- paste("r", 1:length(data.count.SNP.filt$Chromosome), sep="")

# Name each BP 
data.count.SNP.filt$BP <- data.count.SNP.filt$Position

# Get pvalues
data.count.SNP.filt$P <- pchisq(data.count.SNP.filt$LRT, df=1, lower=F)

# Filter data (none of these in my data)
#data.binary.filt <- data.binary.filt[-c(which(data.binary.filt$pvalue == "NaN" ),
#                                        which(data.binary.filt$pvalue == "Inf"),
#                                        which(data.binary.filt$LRT == "Inf")),]

# ================================================================================== #

# Make manhattan plot
manhattan(data.count.SNP.filt, chr="CHR", bp="Position", p="P", 
          ylim = c(0, 350), cex = 0.6, suggestiveline = F, genomewideline = F, xlab = "Position")

# Look at qq-plot of pvalues to check model fit
qqnorm(data.count.SNP.filt$P)

# ================================================================================== #






