## Plotting with plot_admix_function.R code

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
library(qqman)
library(ggplot2)
library(data.table)
library(tidyverse)
library(foreach)

# ================================================================================== #

# Load data 
data.binary.SNP <- read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2.binary.SNPs.gwas.lrt0", header = T, sep = "\t")
str(data.binary.SNP)

# Create unique Chromosome number 
Chr.unique <- unique(data.binary.SNP$Chromosome)
data.binary.SNP$CHR <- as.numeric(factor(data.binary.SNP$Chromosome, levels = Chr.unique))

# Clean data
# Remove LRT values that are -999 (i.e., Sites that fails one of the filters) and are negative
data.binary.SNP.filt <- data.binary.SNP[-c(which(data.binary.SNP$LRT == -999), which(data.binary.SNP$LRT <= 0)), ]

hist(data.binary.SNP.filt$LRT, breaks = 50)

# ================================================================================== #

# Prepare data

# Name each SNP 
data.binary.SNP.filt$SNP <- paste("r", 1:length(data.binary.SNP.filt$Chromosome), sep="")

# Name each BP 
data.binary.SNP.filt$BP <- data.binary.SNP.filt$Position

# Get pvalues
data.binary.SNP.filt$P <- pchisq(data.binary.SNP.filt$LRT, df=1, lower=F)

# ================================================================================== #

# Make manhattan plot
manhattan(data.binary.SNP.filt, chr="CHR", bp="Position", p="P")

# Look at qq-plot of pvalues to check model fit
qqnorm(data.binary.SNP.filt$P)

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

#axisdf = data.binary.filt.graph %>%
#  group_by(CHR) %>%
#  summarize(center=(max(Positioncum) + min(Positioncum) ) / 2 )

ggplot(data.binary.SNP.filt.graph, aes(x=Position.cummulative, y=-log10(P))) +
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

# ================================================================================== #

# Rank-normalize p-values

data.binary.SNP.filt.rn <- data.binary.SNP.filt
data.binary.SNP.filt.rn$rank <- rank(data.binary.SNP.filt$P)
Lp <- length(data.binary.SNP.filt.rn$P)
data.binary.SNP.filt.rn$rn_p_r <- data.binary.SNP.filt.rn$rank/Lp


ggplot(data.binary.SNP.filt.rn, aes(y=-log10(rn_p_r), x=CHR)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  ylab("-log(p)") + xlab("Position") +
  theme_bw()

# ================================================================================== #

# Window analysis

# Define windows
win.bp <- 1e5
step.bp <- 5e4


# Create windows (note: only windows with the number of SNPs in that window >= 30)
wins <- foreach(Chromosome.i=unique(data.binary.SNP.filt.rn$Chromosome),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  
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
              max.p=max(p),
              nSNPs = n(),
              sum.rnp=sum(rn_p_r<=pr.i),
    )  -> win.out
}

# Graph 

win.out %>%
  ggplot(aes(
    x=max.p,
    y=-log10(rnp.binom.p)
  )) + geom_point() + 
  geom_hline(yintercept = -log10(0.05))


win.out %>%
  ggplot(aes(
    x=pos_mean,
    y=-log10(rnp.binom.p)
  )) + geom_point() + geom_hline(yintercept = -log10(0.05)) + theme_bw()



# Create unique Chromosome number 
Chr.unique <- unique(win.out$Chromosome)
win.out$Chr.unique <- as.numeric(factor(win.out$Chromosome, levels = Chr.unique))

ggplot(win.out, aes(y=-log10(rnp.binom.p), x=Chr.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  geom_hline(yintercept = -log10(0.05), color="red") +
  theme_bw()



# ================================================================================== #
