## Plotting with Zoon plot_admix R code

source("~/Desktop/PD_stuffies/scripts/plot_admixture_v5_function.R")

#### K=2
dir=("~/Desktop/PD_stuffies/R_work/admix_K/ALL_pruned/") # path to input files

inName="ALL_2_maf0.05_pctind0.5_maxdepth15_pruned_K2_run2.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run

npops=2

samples <- read.table("~/Desktop/PD_stuffies/R_work/admix_all_pruned_OLD/all.filt.info.txt")

i2p=samples[,1:2] # select the sample name and population name columns
names(i2p)=c("ind","pop") # give column names
row.names(i2p)=i2p$ind # give row names

tbl=read.table(paste(dir,inName,sep=""),header=F) # read in the dataframe from NGSadmix that has % likelihood clustering
row.names(tbl)= i2p$ind # make row names aline with pop dataframe

in.tbl<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL) # merge admix and pop data

# make pop into factor (otherwise will be plotted based in alphabetical order)
in.tbl$pop=factor(in.tbl$pop, levels=c("BA", "BT", "SR", "PB", "CB", "SC", "SH", "IP", "WA", "CA", "GP", "VN", "SB", "CR", "HP", "DB", "BB", "FR", "KR"))

# indicate colors
colors=c("#36648B", "#CD69C9")

# plot
ords=plotAdmixture(data=in.tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)