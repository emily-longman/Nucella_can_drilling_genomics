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
data.binary.SNP <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.binary.all.SNPs.gwas.lrt0", header = T, sep = "\t")
# Pruned SNP list
#data.binary.SNP <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.binary.SNPs.gwas.lrt0", header = T, sep = "\t")
str(data.binary.SNP)

# Create unique Chromosome number 
Chr.unique <- unique(data.binary.SNP$Chromosome)
data.binary.SNP$CHR <- as.numeric(factor(data.binary.SNP$Chromosome, levels = Chr.unique))

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fails one of the filters) and are negative
# Percent of SNPs that fail filter: 48% for full SNP list
dim(data.binary.SNP[-c(which(data.binary.SNP$LRT == -999), which(data.binary.SNP$LRT <= 0)),])[1]/dim(data.binary.SNP)[1]
# Remove sites that fail filter
data.binary.SNP.filt <- data.binary.SNP[-c(which(data.binary.SNP$LRT == -999), which(data.binary.SNP$LRT <= 0)), ]
# Number of sites remaining: 92375 for full SNP list
dim(data.binary.SNP.filt)[1]

hist(data.binary.SNP.filt$LRT, breaks = 50)

# How many SNPs are on each contig:
data.chr.sum <- as.data.frame(table(data.binary.SNP.filt$CHR))

# ================================================================================== #

# Prepare data for graphing using qqman package

# Name each SNP 
data.binary.SNP.filt$SNP <- paste("r", 1:length(data.binary.SNP.filt$Chromosome), sep="")

# Name each BP 
data.binary.SNP.filt$BP <- data.binary.SNP.filt$Position

# Get pvalues
data.binary.SNP.filt$P <- pchisq(data.binary.SNP.filt$LRT, df=1, lower=F)

# ================================================================================== #

# Graph data using qqman package

# Make manhattan plot
manhattan(data.binary.SNP.filt, chr="CHR", bp="Position", p="P")

# Note: genome-wide significance default:  -log10(5e-8) 
# Note: the suggestive line default: -log10(1e-5)
# Typically the genome-wide significance line corresponds to Bonferonni-corrected p-value namely 0.05 divided by the number of SNPs tested

# Regraph adding my own genome-wide significance line
manhattan(data.binary.SNP.filt, chr="CHR", bp="Position", p="P", 
          ylim=c(0,7), suggestiveline = F, genomewideline = F, xlab="Position")
abline(h=-log10(0.05/dim(data.binary.SNP.filt)[1]), col="red")

# Look at qq-plot of pvalues to check model fit
qq(data.binary.SNP.filt$P)

# ================================================================================== #

# Make manhattan plot prettier (https://r-graph-gallery.com/101_Manhattan_plot.html)

data.binary.SNP.filt.graph <- data.binary.SNP.filt %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(data.binary.SNP.filt, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, Position) %>%
  mutate(Position.cummulative=Position+tot)


ggplot(data.binary.SNP.filt.graph, aes(x=Position.cummulative, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), dim(data.binary.SNP.filt)[1] )) +
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

# ================================================================================== #

# Rank-normalize p-values

data.binary.SNP.filt.rn <- data.binary.SNP.filt
data.binary.SNP.filt.rn$rank <- rank(data.binary.SNP.filt$P)
Lp <- length(data.binary.SNP.filt.rn$P)
data.binary.SNP.filt.rn$rn_p_r <- data.binary.SNP.filt.rn$rank/Lp


ggplot(data.binary.SNP.filt.rn, aes(y=-log10(rn_p_r), x=CHR)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  ylab("-log10(rn_p)") + xlab("Position") +
  theme_bw()

# ================================================================================== #

# Window analysis

# Define windows
win.bp <- 1e5
step.bp <- 5e8

# How many SNPs are on each contig:
ggplot(data.chr.sum, aes(x=Freq)) + geom_density() + xlim(0,150)
# Use this information to determine level to filter for number of SNPs in a given window

# Create windows (note: only windows with the number of SNPs in that window >= 5)
wins <- foreach(Chromosome.i=unique(data.binary.SNP.filt.rn$Chromosome),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(Chromosome.i)
                  
                  tmp <- data.binary.SNP.filt.rn %>%
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
  
  
  win.tmp <- data.binary.SNP.filt.rn %>%
    filter(Chromosome == wins[win.i]$Chromosome) %>%
    filter(Position >= wins[win.i]$start & Position <= wins[win.i]$end)
  
  pr.i <- c(0.01)
  
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
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()

# Graph based on position
ggplot(win.out, aes(y=-log10(rnp.binom.p), x=pos_mean/1e6)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()


# ================================================================================== #


# Identify contigs with highly significant rnp p
win.out.sig <- win.out[which(-log10(win.out$rnp.binom.p)>-log10(0.01)),]


# Graph contigs with two significant windows
ggplot(data.binary.SNP.filt.rn[which(data.binary.SNP.filt.rn$Chromosome=="ntLink_5618"),], 
       aes(x=Position/1e6, y=-log10(P))) + 
  geom_line() + theme_classic() + 
  geom_hline(yintercept=-log10(0.01), color="red", linetype="dashed") + ggtitle("ntLink_5618")

ggplot(data.binary.SNP.filt.rn[which(data.binary.SNP.filt.rn$Chromosome=="ntLink_241"),], 
       aes(x=Position/1e6, y=-log10(P))) + 
  geom_line() + theme_classic() + 
  geom_hline(yintercept=-log10(0.01), color="red", linetype="dashed") + ggtitle("ntLink_241")

ggplot(data.binary.SNP.filt.rn[which(data.binary.SNP.filt.rn$Chromosome=="ntLink_4575"),], 
       aes(x=Position/1e6, y=-log10(P))) + 
  geom_line() + theme_classic() + 
  geom_hline(yintercept=-log10(0.01), color="red", linetype="dashed")

# ================================================================================== #

# Create SNP list 
SNPs.Interest <- foreach(i=1:dim(win.out.sig)[1], .combine = "rbind")%do%{
  tmp.snps <- data.binary.SNP.filt %>%
  filter(Chromosome == win.out.sig[i,]$Chromosome) %>%
  filter(Position >= win.out.sig[i,]$pos_min & Position <= win.out.sig[i,]$pos_max)
}

# ================================================================================== #

# Write file of outlier SNPs
write.csv(SNPs.Interest, "Nucella_GWAS_outlier_SNPs.csv")





