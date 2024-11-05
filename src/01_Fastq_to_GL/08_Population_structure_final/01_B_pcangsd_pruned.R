# PCAngsd

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Note: to run script move the covmatrix to the pcangsd folder in results > stats

# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

# Get metadata
metdata_path_from_root <- find_root_file("results", criterion = has_file("README.md"))
setwd(metdata_path_from_root) # Set working directory as path from root
metadata <- read.csv("Metadata.csv", header=T)

# Set relative path of results directory from root
dir(find_root_file("results", criterion = has_file("README.md")))
rel_path_from_root <- find_root_file("results", "stats", "pcangsd_pruned", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(rel_path_from_root)

# Set working directory as path from root
setwd(rel_path_from_root)

# ================================================================================== #

# Load packages
library(ggplot2)
library(ggpubr)

# Load data - cov matrix based on pruned SNP list
COV <- as.matrix(read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2_pval1e6_pruned.cov")) # Read in the genetic covariance matrix

# Extract the principal components from the COV matrix
PCA <- eigen(COV) 

## How much variance is explained by the first few PCs?
var <- round(PCA$values/sum(PCA$values),3)
var[1:5] # First 5 PCs


# A "screeplot" of the eigenvalues of the PCA:
barplot(var, xlab="Eigenvalues of the PCA", ylab="Proportion of variance explained")

## Bring in the bam.list file and extract the sample info:
names <- read.table("Nucella_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".Lanes_merged.bam"))
split = strsplit(names, "-")
#pops <- data.frame(names[1:95], do.call(rbind, split[1:95]))
#names(pops) = c("Ind", "Pop", "Row", "Col")

# Format metadata
str(metadata)
metadata$Total.Drilled <- as.character(metadata$Total.Drilled)

## Graph PCA using ggplot

data=as.data.frame(PCA$vectors)
data=data[,c(1:10)]
data= cbind(data, metadata)

cols=c("#377eB8","#EE9B00", "#7EA16B")

ggscatter(data, x = "V1", y = "V2",
          color = "Site",
          ellipse = T, ellipse.level = 0.95, ellipse.alpha = 0) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols)) 
ggsave("N.canaliculata_Collection_site_pruned_PC1_PC2.jpeg", width = 8, height = 6, device='jpeg', dpi=300)


cols2=c("#94D2BD", "#6d597a")

ggscatter(data, x = "V1", y = "V2",
          color = "Drilled.Binary",
          ellipse = T, ellipse.level = 0.95, ellipse.alpha = 0) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols2)) 
ggsave("N.canaliculata_Drilled_Binary_pruned_PC1_PC2.jpeg", width = 8, height = 6, device='jpeg', dpi=300)

ggscatter(data, x = "V1", y = "V3",
          color = "Drilled.Binary",
          ellipse = T, ellipse.level = 0.95, ellipse.alpha = 0) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC3: (",var[3]*100,"%)")) +
  scale_color_manual(values=c(cols2)) 
ggsave("N.canaliculata_Drilled_Binary_pruned_PC1_PC3.jpeg", width = 8, height = 6, device='jpeg', dpi=300)


cols3=c("#377eB8","#EE9B00","#94D2BD","#AE2012", "#6d597a", "#7EA16B")

ggscatter(data, x = "V1", y = "V2",
          color = "Total.Drilled",
          ellipse = T, ellipse.level = 0.95, ellipse.alpha = 0) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols3)) 
ggsave("N.canaliculata_Total_Drilled_pruned_PC1_PC2.jpeg", width = 8, height = 6, device='jpeg', dpi=300)
