argv <- commandArgs(T)
infile <- argv[1] 

# Read in bamlist
bam_list<-read.table(infile, header=F)

# Sample 7 of the 11 individuals

bam_list_List <- list()

for(i in 1:100){
    bam_list_List[[i]] <- bam_list[sample(11, 7, replace=F),]
}

for(i in seq_along(bam_list_List)) {
    write.table(bam_list_List[[i]], paste("Olaa_bam_filelist_reduced_", i, ".list", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
}