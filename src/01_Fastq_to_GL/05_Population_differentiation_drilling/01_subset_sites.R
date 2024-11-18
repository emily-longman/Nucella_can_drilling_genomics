# This R script will subset the bamfile lists so that the bam lists for drilling will have the same sample size (i.e., n = 75)

argv <- commandArgs(T)
Drilling_bamlist <- argv[1]
Drilling_group <- argv[2]

# Read in bamlist
Drilling_bamlist<-read.table(Drilling_bamlist, header=F)
n_group<-dim(Drilling_bamlist)[1]

# Need to sample for 75 individuals
n_ind <- 75

bam_list <- Drilling_bamlist[sample(n_group, 75, replace=F),]

# Write subset bam list 
write.table(bam_list, paste("Nucella_bam_", Drilling_group, "_subset.list", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)