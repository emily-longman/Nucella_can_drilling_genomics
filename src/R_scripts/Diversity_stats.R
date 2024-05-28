# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))
# Set relative path from root
rel_path_from_root <- find_root_file("results", "diversity_stats", criterion = has_file("README.md"))

# List files in this folder to make sure you're in the right spot.
list.files(rel_path_from_root)
# Set working directory as path from root
setwd(rel_path_from_root)

# Load data
theta.FB <- read.table("FB_SNPs.thetas", sep="\t", header=T)
theta.HC <- read.table("FB_SNPs.thetas", sep="\t", header=T)
theta.MP <- read.table("FB_SNPs.thetas", sep="\t", header=T)


theta.FB$tWsite = theta.FB$tW/theta.FB$nSites #scales the theta-W by the number of sites
theta.FB$tPsite = theta.FB$tP/theta.FB$nSites #scales the theta-Pi by the number of sites

summary(theta.FB)

# You can order the contig list to show you the contigs with the highest values of Tajima's D, or the lowest

head(theta.FB[order(theta.FB$Tajima, decreasing = TRUE),]) # top 10 Tajima's D values

head(theta.FB[order(theta.FB$Tajima, decreasing = FALSE),]) # bottom 10 Tajima's D values

#You can also look for contigs that have combinations of high Tajima's D and low diversity -- these could represent outliers for selection
#theta[which(theta$Tajima>1.5 & theta$tPsite<0.001),]

sfs.FB<-scan('FB_SNPs.sfs')
sfs.FB<-sfs.FB[-c(1,which(sfs.FB==0))]
sfs.FB<-sfs.FB/sum(sfs.FB)

# Be sure to replace "9999" with your pop code in the "main" legend below
barplot(sfs.FB, xlab="Contigs",
        names=1:length(sfs.FB),
        ylab="Proportions",
        main="Pop FB Site Frequency Spectrum",
        col='blue')

# Put the nucleotide diversities, Tajima's D, and SFS into a 4-panel figure
par(mfrow=c(2,2))
hist(theta.FB$tWsite, xlab="theta-W", main="Watterson's theta")
hist(theta.FB$tPsite, xlab="theta-Pi", main="Pairwise Nucleotide Diversity")
hist(theta.FB$Tajima, xlab="D", main="Tajima's D")
barplot(sfs.FB,names=1:length(sfs.FB),main='Site Frequency Spectrum')


# To reset the panel plotting, execute the line below:
dev.off()