# Plot GWAS results

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
library(ggplot2)
library(dplyr)

# ================================================================================== #

# UPDATE WITH NEW lrt file - this one included regions
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
manhattan(data.binary.filt, chr="CHR", bp="Position", p="P")

# Look at qq-plot of pvalues to check model fit
qqnorm(data.binary.filt$P)

# ================================================================================== #

# Highlight values that exceed a threshold - i.e., the highest LRT value from a random phenotype test 

# Load data
data.binary.random <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.binary.RANDOM.gwas.lrt0", header = T, sep = "\t")
str(data.binary.random)

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fails one of the filters) and are negative
data.binary.random.filt <- data.binary.random[-c(which(data.binary.random$LRT == -999), which(data.binary.random$LRT <= 0)), ]

# Summarize LRT of random data
summary(data.binary.random.filt$LRT, na.rm = T)
max(data.binary.random.filt$LRT, na.rm = T)

# ================================================================================== #

# Make a list of the candidate loci
candidates <- data.binary.filt[which(data.binary.filt$LRT > 20),]$SNP

# Graph manhattan and highlight candidate loci
manhattan(data.binary.filt, chr="CHR", bp="Position", p="P", highlight=candidates, cex=0.6, col="grey")
manhattan(data.binary.filt, chr="CHR", bp="Position", p="P", annotatePval = 0.0001)

# ================================================================================== #

# Make manhattan plot prettier (https://r-graph-gallery.com/101_Manhattan_plot.html)

data.binary.filt.graph <- data.binary.filt %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(data.binary.filt, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, Position) %>%
  mutate(Position.cummulative=Position+tot)

#axisdf = data.binary.filt.graph %>%
#  group_by(CHR) %>%
#  summarize(center=(max(Positioncum) + min(Positioncum) ) / 2 )

ggplot(data.binary.filt.graph, aes(x=Position.cummulative, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 5534 )) +
  # custom X axis:
  #scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(limits=c(0, 6), expand = c(0, 0) ) + # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )



