# PCAngsd

# ================================================================================== #

# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

dir(find_root_file("results", criterion = has_file("README.md")))
# Set relative path from root
rel_path_from_root <- find_root_file("results", "stats", "pcangsd", criterion = has_file("README.md"))

# List files in this folder to make sure you're in the right spot.
list.files(rel_path_from_root)
# Set working directory as path from root
setwd(rel_path_from_root)

# ================================================================================== #

# Clear memory
rm(list=ls()) 

# Load packages
library(ggplot2)
library(ggpubr)

# Load data
COV <- as.matrix(read.table("Nucella_poly_covmatrix.cov")) # Read in the genetic covariance matrix

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

# Load metadata

metadata <- read.csv("Metadata.csv", header=T)
str(metadata)
metadata$Total.Drilled <- as.character(metadata$Total.Drilled)

## A more beautiful PCA plot using ggplot :)

data=as.data.frame(PCA$vectors)
data=data[,c(1:10)]
data= cbind(data, metadata)

cols=c("#377eB8","#EE9B00", "#7EA16B")
#cols=c("#377eB8","#EE9B00","#0A9396","#94D2BD","#FFCB69","#005f73","#E26D5C","#AE2012", "#6d597a", "#7EA16B","#d4e09b", "gray70")

ggscatter(data, x = "V1", y = "V2",
          color = "Site",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))



cols2=c("#94D2BD", "#6d597a")

ggscatter(data, x = "V1", y = "V2",
          color = "Drilled.Binary",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols2), name="Source population") +
  guides(colour = guide_legend(nrow = 2))

ggscatter(data, x = "V1", y = "V3",
          color = "Drilled.Binary",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols2), name="Source population") +
  guides(colour = guide_legend(nrow = 2))


cols3=c("#377eB8","#EE9B00","#94D2BD","#AE2012", "#6d597a", "#7EA16B")

ggscatter(data, x = "V1", y = "V3",
          color = "Total.Drilled",
          mean.point = TRUE,
          star.plot = TRUE) +
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols3), name="Source population") +
  guides(colour = guide_legend(nrow = 2))
