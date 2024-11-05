# This R script will subset the bamfile lists so that the sites have the same sample size (i.e., n = 59)

argv <- commandArgs(T)
site_name <- argv[1]
site <- argv[2]

# Read in bamlist
site_bamlist<-read.table(site, header=F)
n_site<-dim(site_bamlist)[1]

# Need to sample for 59 individuals
n_ind <- 59

bam_list <- site_bamlist[sample(n_site, 59, replace=F),]

# Write subset bam list 
write.table(bam_list, paste(site_name, "_bam_subset", ".list", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)