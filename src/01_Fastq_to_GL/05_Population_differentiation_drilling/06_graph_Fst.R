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
library(data.table)
library(foreach)

# ================================================================================== #

# Load data 
data <- read.table("Drilled_Not.Drilled_maf0.05_pctind0.5_mindepth0.3_maxdepth2_subset_nMAF.slidingwindow", header = T, row.names=NULL)
# Fix column names
colnames(data) <- c("region", "chr", "midPos", "Nsites", "Fst")

# Create unique Chromosome number 
Chr.unique <- unique(data$chr)
data$CHR.unique <- as.numeric(factor(data$chr, levels = Chr.unique))

# Graph Fst against chromosome 
ggplot(data, aes(y=Fst, x=CHR.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  ylab("Window Fst") + xlab("Position") +
  theme_bw()

# ================================================================================== #

# Rank-normalize Fst values

data.rn <- data
data.rn$rank <- rank(-data$Fst)
Lfst <- length(data$Fst)
data.rn$rn_fst_r <- data.rn$rank/Lfst

ggplot(data.rn, aes(y=-log10(rn_fst_r), x=CHR.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  ylab("Window Fst") + xlab("Position") +
  theme_bw()

# Window analysis

win.bp <- 25000
step.bp <- 5000


wins <- foreach(chr.i=unique(data.rn$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  
                  tmp <- data.rn %>%
                    filter(chr == chr.i)
                  
                  S=dim(tmp)[1]
                
                  if(S >= 30){
                    o =
                    data.table(chr=chr.i,
                               S=dim(tmp)[1],
                               start=seq(from=min(tmp$midPos), to=max(tmp$midPos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$midPos), to=max(tmp$midPos)-win.bp, by=step.bp) + win.bp)
                    return(o)
                    
                  }   
                  else {message("fails S filter")}
                }

wins[,i:=1:dim(wins)[1]]
####
dim(wins)


####


### start the summarization process
win.out <- foreach(win.i=1:dim(wins)[1], 
#win.out <- foreach(win.i=1:10000, 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
  win.tmp <- data.rn %>%
    filter(chr == wins[win.i]$chr) %>%
    filter(midPos >= wins[win.i]$start & midPos <= wins[win.i]$end)
  
  pr.i <- c(0.05)
  
  win.tmp %>% 
    filter(!is.na(rn_fst_r)) %>%
    summarise(chr = wins[win.i]$chr,
              pos_mean = mean(midPos),
              pos_mean = mean(midPos),
              pos_min = min(midPos),
              pos_max = max(midPos),
              win=win.i,
              pr=pr.i,
              rnp.pr=c(mean(rn_fst_r<=pr.i)),
              rnp.binom.p=c(binom.test(sum(rn_fst_r<=pr.i), 
                                       length(rn_fst_r), pr.i)$p.value),
              max.fst=max(Fst),
              nSNPs = n(),
              sum.rnp=sum(rn_fst_r<=pr.i),
    )  -> win.out
  }



win.out %>%
  ggplot(aes(
    x=max.fst,
    y=-log10(rnp.binom.p)
  )) + geom_point() + 
  geom_hline(yintercept = -log10(0.05))
