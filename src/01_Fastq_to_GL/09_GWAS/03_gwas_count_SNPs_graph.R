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
library(data.table)
library(tidyverse)
library(foreach)

# ================================================================================== #

# Load data 
# All SNPs
data.count.SNP <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.count.all.SNPs.gwas.lrt0", header = T, sep = "\t")
# Pruned SNP list
#data.count.SNP <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.count.SNPs.gwas.lrt0", header = T, sep = "\t")
str(data.count.SNP)

# Create unique Chromosome number 
Chr.unique <- unique(data.count.SNP$Chromosome)
data.count.SNP$CHR <- as.numeric(factor(data.count.SNP$Chromosome, levels = Chr.unique))

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fail one of the filters) and are negative
# Percent of SNPs that fail filter: 6.7% for full SNP list
dim(data.count.SNP[-c(which(data.count.SNP$LRT == -999), which(data.count.SNP$LRT <= 0)),])[1]/dim(data.count.SNP)[1]
# Remove sites that fail filter
data.count.SNP.filt <- data.count.SNP[-c(which(data.count.SNP$LRT == -999), which(data.count.SNP$LRT <= 0)), ]
# Number of sites remaining: 92375 for full SNP list
dim(data.count.SNP.filt)[1]

hist(data.count.SNP.filt$LRT, breaks = 50)

# How many SNPs are on each contig:
data.chr.sum <- as.data.frame(table(data.count.SNP.filt$CHR))

# ================================================================================== #

# Prepare data

# Name each SNP 
data.count.SNP.filt$SNP <- paste("r", 1:length(data.count.SNP.filt$Chromosome), sep="")

# Name each BP 
data.count.SNP.filt$BP <- data.count.SNP.filt$Position

# Get pvalues
data.count.SNP.filt$P <- pchisq(data.count.SNP.filt$LRT, df=1, lower=F)

# Filter data (none in this dataset)
#data.count.SNP.filt.tmp <- data.count.SNP.filt[-c(which(data.count.SNP.filt$P == "NaN"),
#                                        which(data.count.SNP.filt$P == "Inf"),
#                                        which(data.count.SNP.filt$LRT == "Inf")), ]

# ================================================================================== #

# Make manhattan plot
manhattan(data.count.SNP.filt, chr="CHR", bp="Position", p="P", ylim = c(0, 400), xlab = "Position")

# Note: genome-wide significance default:  -log10(5e-8) 
# Note: the suggestive line default: -log10(1e-5)
# Typically the genome-wide significance line corresponds to Bonferonni-corrected p-value namely 0.05 divided by the number of SNPs tested

# Regraph adding my own genome-wide significance line
manhattan(data.count.SNP.filt, chr="CHR", bp="Position", p="P", 
          ylim=c(0,400), suggestiveline = F, genomewideline = F, xlab="Position")
abline(h=-log10(0.05/dim(data.count.SNP.filt)[1]), col="red")

# Look at qq-plot of pvalues to check model fit
qq(data.count.SNP.filt$P)

# ================================================================================== #

# NOTE: Problem - getting some SNPs with p = 0

# Make manhattan plot prettier (https://r-graph-gallery.com/101_Manhattan_plot.html)

data.count.SNP.filt.graph <- data.count.SNP.filt %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(data.count.SNP.filt, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, Position) %>%
  mutate(Position.cummulative=Position+tot)

#axisdf = data.binary.filt.graph %>%
#  group_by(CHR) %>%
#  summarize(center=(max(Positioncum) + min(Positioncum) ) / 2 )

ggplot(data.count.SNP.filt.graph, aes(x=Position.cummulative, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 5534)) +
  # custom X axis:
  #scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(limits=c(0, 400), expand = c(0, 0) ) + # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# ================================================================================== #

# Rank-normalize p-values

data.count.SNP.filt.rn <- data.count.SNP.filt
data.count.SNP.filt.rn$rank <- rank(data.count.SNP.filt$P)
Lp <- length(data.count.SNP.filt.rn$P)
data.count.SNP.filt.rn$rn_p_r <- data.count.SNP.filt.rn$rank/Lp


ggplot(data.count.SNP.filt.rn, aes(y=-log10(rn_p_r), x=CHR)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  ylab("-log(p)") + xlab("Position") +
  theme_bw()

# ================================================================================== #

# Window analysis

# Define windows
win.bp <- 1e5
step.bp <- 5e4

# How many SNPs are on each contig:
ggplot(data.chr.sum, aes(x=Freq)) + geom_density() + xlim(0,150)
# Use this information to determine level to filter for number of SNPs in a given window

# Create windows (note: only windows with the number of SNPs in that window >= 5)
wins <- foreach(Chromosome.i=unique(data.count.SNP.filt.rn$Chromosome),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(Chromosome.i)
                  
                  tmp <- data.count.SNP.filt.rn %>%
                    filter(Chromosome == Chromosome.i)
                  
                  S=dim(tmp)[1]
                  
                  if(S >= 5){
                    o =
                      data.table(Chromosome=Chromosome.i,
                                 S=dim(tmp)[1],
                                 start=seq(from=min(tmp$Position), to=max(tmp$Position)-win.bp, by=step.bp),
                                 end=seq(from=min(tmp$Position), to=max(tmp$Position)-win.bp, by=step.bp) + win.bp)
                    return(o)
                    
                  }   
                  else {message("fails S filter")}
                }

wins[,i:=1:dim(wins)[1]]
dim(wins)

# ================================================================================== #

# Start the summarization process
win.out <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
  win.tmp <- data.count.SNP.filt.rn %>%
    filter(Chromosome == wins[win.i]$Chromosome) %>%
    filter(Position >= wins[win.i]$start & Position <= wins[win.i]$end)
  
  pr.i <- c(0.05)
  
  win.tmp %>% 
    filter(!is.na(rn_p_r)) %>%
    summarise(Chromosome = wins[win.i]$Chromosome,
              pos_mean = mean(Position),
              pos_mean = mean(Position),
              pos_min = min(Position),
              pos_max = max(Position),
              win=win.i,
              pr=pr.i,
              rnp.pr=c(mean(rn_p_r<=pr.i)),
              rnp.binom.p=c(binom.test(sum(rn_p_r<=pr.i), 
                                       length(rn_p_r), pr.i)$p.value),
              max.p=max(P),
              nSNPs = n(),
              sum.rnp=sum(rn_p_r<=pr.i),
    )  -> win.out
}

# ================================================================================== #

# Graph 

# Create unique Chromosome number 
Chr.unique <- unique(win.out$Chromosome)
win.out$Chr.unique <- as.numeric(factor(win.out$Chromosome, levels = Chr.unique))

ggplot(win.out, aes(y=-log10(rnp.binom.p), x=Chr.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  geom_hline(yintercept = -log10(0.05), color="red") +
  theme_bw()


