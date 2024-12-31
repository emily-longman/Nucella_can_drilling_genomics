# Plotting with plot_admix_function.R code

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Note: to run script move the ngs_admix folder in results 

# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

# Set relative path of results directory from root
dir(find_root_file("results", criterion = has_file("README.md")))
results_path_from_root <- find_root_file("results", "stats", "fst", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(results_path_from_root)

# Set working directory as path from root
setwd(results_path_from_root)

# ================================================================================== #

# Graph fst along sliding window

# ================================================================================== #

library(ggplot2)
library(data.table)
library(foreach)

# ================================================================================== #

# Load data 
data <- read.table("Drilled_Not.Drilled_maf0.05_pctind0.5_mindepth0.3_maxdepth2_subset_nMAF.slidingwindow", header = T, row.names=NULL)
# Fix column names
colnames(data) <- c("region", "chr", "midPos", "Nsites", "Fst")

# Create unique Chromosome number 
Chr.unique <- unique(data$chr)
data$CHR.unique <- as.numeric(factor(data$chr, levels = Chr.unique))

# Graph Fst against chromosome 
ggplot(data, aes(y=Fst, x=CHR.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  ylab("Window Fst") + xlab("Position") +
  theme_bw()

# ================================================================================== #

# Rank-normalize Fst values

data.rn <- data
data.rn$rn_fst <- rank(data$Fst)/length(data$Fst)
data.rn$rn_fst_r <- rank(data$Fst, ties="random")/length(data$Fst)

# Window analysis

win.bp <- 25000
step.bp <- 5000


