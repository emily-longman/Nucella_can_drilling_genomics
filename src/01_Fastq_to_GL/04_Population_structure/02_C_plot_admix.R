## Plotting with plot_admix_function.R code

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Note: to run script move the ngs_admix folder in results 

# Set path as main Github repo
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

# Get metadata
metdata_path_from_root <- find_root_file("results", criterion = has_file("README.md"))
setwd(metdata_path_from_root) # Set working directory as path from root
samples <- read.csv("Metadata.csv", header=T)

# Set relative path of results directory from root
dir(find_root_file("results", criterion = has_file("README.md")))
results_path_from_root <- find_root_file("results", "stats", "ngs_admix", "K_output", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(results_path_from_root)

# Set working directory as path from root
setwd(results_path_from_root)

# ================================================================================== #

# Look at log file produced in step "02_B_admix_all". Identify which K has lowest likelihood.  
# Further, if curious, plot the iteration for a given K with the lowest likelihood

# ================================================================================== #

# Load packages
library(stringr)
library(ggplot2)
library(colorspace)
library(RColorBrewer)

# ================================================================================== #

# Get R plot admix function
source("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/01_Fastq_to_GL/04_Population_structure/02_C_plot_admix_function.R")

# ================================================================================== #

# Choose K using two methods - maximum likelihood and CLUMPAK

# ================================================================================== #

# Identify the list of log files
data <- list.files(pattern = ".log", full.names = T)

# Use lapply to read in all our log files at once
bigData<-lapply(1:30, FUN = function(i) readLines(data[i]))

# Extract the line that starts with "b" from each file and return it as a list
foundset<-sapply(1:30, FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])

# Extract the first number in the string. Do this with the function sub
as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )

# Make a dataframe with an index 1:3, this corresponds to our K values
logs <- data.frame(K = rep(1:3, each=10))

# Add to it our likelihood values
logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

# Calculate our delta K and probability
tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

# --> choose K value with highest maximum likelihood

# ================================================================================== #

# Get K with CLUMPAK (https://clumpak.tau.ac.il/bestK.html)
# However, since all K=1 have the same likelihood can't use CLUMPAK 
# "Cannot calculate best k because Standard deviation of 'Ln Prob of Data' for K=1 is zero."

# Change path to load in log file
results_path_from_root <- find_root_file("results", "stats", "ngs_admix", criterion = has_file("README.md"))
setwd(results_path_from_root)

# Import logfile & convert to CLUMPAK format
logs <- as.data.frame(read.table("logfile"))

# Add K
logs$K <- c(rep("1", 10), rep("2", 10), rep("3", 10))

# Write table
write.table(logs[, c(2, 1)], "logfile_formatted", row.names = F, col.names = F, quote = F)

# --> The re-formatted logfile can be used with Clumpak to determine your most likely K.
# Upload the file - check log probability table file

# ================================================================================== #

# Change path for subsequent graphing

# Set relative path of results directory from root
results_path_from_root <- find_root_file("results", "stats", "ngs_admix", "K_output", criterion = has_file("README.md"))
setwd(results_path_from_root)

# ================================================================================== #

# Graph using two methods - barplots and plotAdmixture code
# For visualization purposes graph K=2 and K=3

# ================================================================================== #

### Visualize NGSadmix output with structure plots

# K=2

# Read in a qopt file. It is formatted such that there is one row for each individual and one column for the percentage of admixture for each of the K's
q<-read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2_pval1e6_K{2}_run{3}.qopt")

# Indicate colors
colors=sequential_hcl(2, palette = "TealGrn") # Number of colors you want and then pallete style

barplot(t(q),
        col=colors, #change colors to number of K
        names=samples$Site,
        las=2,
        space=0,
        border=NA,
        xlab="Individuals",
        ylab="Admixture proportions for K=2 pruned")
abline(v = c(59, 127), lty = 5, lwd = 2, col = "black")
ggsave("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/stats/ngs_admix/N.canaliculata_Admixture_K=2_Site_barplot.jpeg", width = 8, height = 6, device='jpeg', dpi=300)

####### 

# K=3

# Read in a qopt file. It is formatted such that there is one row for each individual and one column for the percentage of admixture for each of the K's
q<-read.table("Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2_pval1e6_K{3}_run{4}.qopt")

# Indicate colors
colors=sequential_hcl(3, palette = "TealGrn") # Number of colors you want and then pallete style

barplot(t(q),
        col=colors, 
        names=samples$Site,
        las=2,
        space=0,
        border=NA,
        xlab="Individuals",
        ylab="Admixture proportions for K=2 pruned")
abline(v = c(59, 127), lty = 5, lwd = 2, col = "white")
ggsave("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/stats/ngs_admix/N.canaliculata_Admixture_K=3_Site_barplot.jpeg", width = 8, height = 6, device='jpeg', dpi=300)

# ================================================================================== #

# Graph using plotAdmixture script

#### Set K=2
npops=2 #npops (i,e., number of groups) = K

# Look at log file to find which has lowest 
inName="Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2_pval1e6_K{2}_run{3}.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run

########### 

# Look at Binary Drilling

# Read in table
tbl=read.table(inName, sep="", header = F) # read in the dataframe from NGSadmix that has % likelihood clustering

# Get metadata
i2p=samples[,c(1,4)] # select the sample name and population name columns
names(i2p)=c("ind","drilled") # give column names
row.names(i2p)=i2p$ind # give row names

# Make row names aline with pop dataframe
row.names(tbl)= i2p$ind 

# Merge admix and pop data
in.tbl<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL) 

# Make pop into factor (otherwise will be plotted based in alphabetical order)
in.tbl$pop=factor(in.tbl$drilled, levels=c("Drilled","Not.Drilled"))

# Indicate colors
colors=sequential_hcl(3, palette = "TealGrn") # Number of colors you want and then pallete style

# Plot
ords=plotAdmixture(data=in.tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)
ggsave("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/stats/ngs_admix/N.canaliculata_Admixture_K=2_Drilled.jpeg", width = 8, height = 6, device='jpeg', dpi=300)

########### 

# Look at Collection Source

# Read in table
tbl=read.table(inName, sep="", header = F) # read in the dataframe from NGSadmix that has % likelihood clustering

# Get metadata
i2p=samples[,c(1,2)] # select the sample name and population name columns
names(i2p)=c("ind","source") # give column names
row.names(i2p)=i2p$ind # give row names

# Make row names aline with pop dataframe
row.names(tbl)= i2p$ind 

# Merge admix and pop data
in.tbl<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL) 

# Make source into factor (otherwise will be plotted based in alphabetical order)
in.tbl$pop=factor(in.tbl$source, levels=c("FB", "HC", "MP"))

# Indicate colors
colors=sequential_hcl(3, palette = "TealGrn") # Number of colors you want and then pallete style

# Plot
ords=plotAdmixture(data=in.tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)
ggsave("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/stats/ngs_admix/N.canaliculata_Admixture_K=2_Site.jpeg", width = 8, height = 6, device='jpeg', dpi=300)

# ================================================================================== #


#### Set K=3
npops=3 #npops (i,e., number of groups) = K

# Look at log file to find which has lowest 
inName="Nucella_SNPs_maf0.05_pctind0.5_mindepth0.3_maxdepth2_pval1e6_K{3}_run{4}.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run

########### 

# Look at Binary Drilling

# Read in table
tbl=read.table(inName, sep="", header = F) # read in the dataframe from NGSadmix that has % likelihood clustering

# Get metadata
i2p=samples[,c(1,4)] # select the sample name and population name columns
names(i2p)=c("ind","drilled") # give column names
row.names(i2p)=i2p$ind # give row names

# Make row names aline with pop dataframe
row.names(tbl)= i2p$ind 

# Merge admix and pop data
in.tbl<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL) 

# Make pop into factor (otherwise will be plotted based in alphabetical order)
in.tbl$pop=factor(in.tbl$drilled, levels=c("Drilled","Not.Drilled"))

# Indicate colors
colors=sequential_hcl(3, palette = "TealGrn") # Number of colors you want and then pallete style

# Plot
ords=plotAdmixture(data=in.tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)
ggsave("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/stats/ngs_admix/N.canaliculata_Admixture_K=3_Drilled.jpeg", width = 8, height = 6, device='jpeg', dpi=300)

########### 

# Look at Collection Source

# Read in table
tbl=read.table(inName, sep="", header = F) # read in the dataframe from NGSadmix that has % likelihood clustering

# Get metadata
i2p=samples[,c(1,2)] # select the sample name and population name columns
names(i2p)=c("ind","source") # give column names
row.names(i2p)=i2p$ind # give row names

# Make row names aline with pop dataframe
row.names(tbl)= i2p$ind 

# Merge admix and pop data
in.tbl<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL) 

# Make source into factor (otherwise will be plotted based in alphabetical order)
in.tbl$pop=factor(in.tbl$source, levels=c("FB", "HC", "MP"))

# Indicate colors
colors=sequential_hcl(3, palette = "TealGrn") # Number of colors you want and then pallete style

# Plot
ords=plotAdmixture(data=in.tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)
ggsave("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/stats/ngs_admix/N.canaliculata_Admixture_K=3_Site.jpeg", width = 8, height = 6, device='jpeg', dpi=300)
