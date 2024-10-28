## Plotting with Zoon plot_admix R code

source("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/src/01_Fastq_to_GL/04_Population_structure/04_C_plot_admix_function.R")

#### K=2
dir=("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/ngs_admix/K_output/") # path to input files

inName="Nucella_all_maf_K{4}_run{3}.qopt" # name of the input file to plot, output of ngsAdmix or ADMIXTURE run

npops=4

samples <- read.csv("/Users/emilylongman/Documents/GitHub/Nucella_can_drilling_genomics/results/ngs_admix/Metadata.csv", header=T)

i2p=samples[,c(1,4)] # select the sample name and population name columns
names(i2p)=c("ind","pop") # give column names
row.names(i2p)=i2p$ind # give row names

tbl=read.table(paste(dir,inName,sep=""),header=F) # read in the dataframe from NGSadmix that has % likelihood clustering
row.names(tbl)= i2p$ind # make row names aline with pop dataframe

in.tbl<-transform(merge(tbl,i2p,by=0), row.names=Row.names, Row.names=NULL) # merge admix and pop data

# make pop into factor (otherwise will be plotted based in alphabetical order)
#in.tbl$pop=factor(in.tbl$pop, levels=c("FB", "HC", "MP"))
in.tbl$pop=factor(in.tbl$pop, levels=c("Drilled","Not.Drilled"))

# indicate colors
colors=c("tomato", "lightblue", "wheat","olivedrab")
colors2=c("#36648B", "#CD69C9")

# plot
ords=plotAdmixture(data=in.tbl,npops=npops,grouping.method="distance",vshift=0.1, colors = colors)
