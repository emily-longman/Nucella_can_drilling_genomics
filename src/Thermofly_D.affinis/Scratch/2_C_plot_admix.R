## Plotting with Zoon plot_admix R code

install.packages(c('data.table', 'ggplot2'))
library(data.table)
library(ggplot2)

source("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/Thermofly/4_Pop_structure/2_C_plot_admix_function.R")

#### K=2
dir=("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/Thermofly/K_output/") # path to input files

inName="Thermofly_GL_reduced_minInd_16_depth_4_maf_K{2}_run{2}.fopt.gz" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run

npops=2

samples <- fread("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/Thermofly/Thermofly_metadata.tsv", header=T)
samples_bas<-meta_data[which(samples$species_initial=="bas"),]
# Remove metadata for individuals with low coverage
samples_bas_reduced <- samples_bas[-which(samples_bas$sampleId =="D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_9.1"| 
                                            samples_bas$sampleId == "D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_12.1" |
                                            samples_bas$sampleId == "D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_13.1" |
                                            samples_bas$sampleId == "D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_15.1")]


i2p=samples_bas_reduced[,c(1,15)] # select the sample name and population name columns
names(i2p)=c("ind","pop") # give column names
row.names(i2p)=i2p$ind # give row names

# BELOW IS NOT WORKING!

tbl=read.table(paste(dir,inName,sep=""),header=F) # read in the dataframe from NGSadmix that has % likelihood clustering
row.names(tbl)= i2p$ind # make row names aline with pop dataframe

in.tbl<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL) # merge admix and pop data

# make pop into factor (otherwise will be plotted based in alphabetical order)
#in.tbl$pop=factor(in.tbl$pop, levels=c("FB", "HC", "MP"))
in.tbl$pop=factor(in.tbl$pop, levels=c("Tom's Trail","Olaa Forest"))

# indicate colors
colors=c("tomato", "lightblue", "wheat","olivedrab")
colors2=c("#36648B", "#CD69C9")

# plot
ords=plotAdmixture(data=in.tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)
