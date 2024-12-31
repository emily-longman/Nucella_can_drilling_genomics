# Run the cmh test for allele differences between drilled versus not drilled groups

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
results_path_from_root <- find_root_file("results", "stats", "cmh", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(results_path_from_root)

# Set working directory as path from root
setwd(results_path_from_root)

# ================================================================================== #

# Perform the Cochran-Mantel-Haenszel (CMH) test.

# ================================================================================== #

# https://github.com/ThomasTaus/poolSeq?tab=readme-ov-file
install.packages("/Users/emilylongman/Documents/Software/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")
library(poolSeq)
library(qvalue)

# ================================================================================== #

# Read baypass input file
bpass.in<- read.table("by_group_0.05_pctind0.25_mindepth0.3_maxdepth2.mafs.pruned.baypass", header=F)

# Transform (rows=pops, cols=inds)
cmh.in <- t(bpass.in)
dim(cmh.in)

#get minor allele counts for drilled snails
p0<- cmh.in[c(1),]

#get major allele counts for drilled snails
P0<- cmh.in[c(2),]

#get minor allele counts for non-drilled snails
np0<- cmh.in[c(3),]

#get major allele counts for non-drilled snails
NP0<- cmh.in[c(4),]

# perform CMH-test for all empirical loci
p.values <- cmh.test(A0=P0, a0=p0, At=NP0, at=np0, min.cov=1, max.cov=0.99, min.cnt=1, log=TRUE)

# plot -log10 transformed p-values
plot(p.values, main="CMH Manhattan plot", xlab="SNP", ylab="-log10(p)", pch=".")

all_qvals <- qvalue(p.values)
FDR <- 0.1
which(all_qvals$qvalues < 0.1)
[1] 342056

pval_cutoff <-max(all_qvals$pvalues[all_qvals$qvalues <= FDR])
1.050089e-07 