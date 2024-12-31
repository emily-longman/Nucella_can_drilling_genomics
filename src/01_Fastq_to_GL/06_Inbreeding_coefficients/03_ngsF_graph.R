# Graph inbreeding coefficient across collection locations.

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
results_path_from_root <- find_root_file("results", "stats", "ngsF", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(results_path_from_root)

# Set working directory as path from root
setwd(results_path_from_root)

# ================================================================================== #

# Graph inbreeding coefficient.

# ================================================================================== #

library(ggplot2)
library(tidyverse)
library(colorspace)
library(RColorBrewer)

# ================================================================================== #

# Load data 
FB <- read.table("FB.indF", header = F, row.names=NULL)
HC <- read.table("HC.indF", header = F, row.names=NULL)
MP <- read.table("MP.indF", header = F, row.names=NULL)

# Add column with site name 
FB$Site <- "FB"
HC$Site <- "HC"
MP$Site <- "MP"

# Join datasets 
inbreeding <- dplyr::bind_rows(FB, HC, MP)

# Change column names
colnames(inbreeding) <- c("F", "Site")

# ================================================================================== #

# Boxplot of inbreeding coefficient 

# Indicate colors
cols=c("#377eB8","#EE9B00", "#7EA16B")

# Graph 
ggplot(inbreeding, aes(Site, F, fill = Site)) + geom_boxplot() + scale_fill_manual(values = cols) + 
  geom_jitter(width = 0.2) +
  labs(x="Collection Site", 
       y=expression('Inbreeding Coefficient'~italic((F))), fill="Var1") + 
  theme_classic() + guides(fill="none")
  
# ================================================================================== #

# Model data with a linear model
inbreeding.mod <- lm(F ~ Site, inbreeding)

# Check model assumptions
par(mfrow=c(1,3))
plot(inbreeding.mod, which=2:3)
hist(inbreeding.mod$residuals)

# Model results
anova(inbreeding.mod)

