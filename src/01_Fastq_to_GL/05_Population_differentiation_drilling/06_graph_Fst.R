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

# ================================================================================== #

# Load data 
data <- read.table("Drilled_Not.Drilled_maf0.05_pctind0.5_mindepth0.3_maxdepth2_subset_nMAF.slidingwindow", header = T, row.names=NULL)


# Create unique Chromosome number 
Chr.unique <- unique(data$region)
data$CHR.unique <- as.numeric(factor(data$region, levels = Chr.unique))

ggplot(data, aes(y=Nsites, x=CHR.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  theme_bw() 
