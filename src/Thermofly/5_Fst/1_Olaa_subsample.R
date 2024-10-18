argv <- commandArgs(T)
infile <- argv[1] 

# Read in bamlist
bam_list<-read.table(infile, header=F)

# Sample 7 of the 11 individuals
bam_list_i <- matrix(ncol=1, nrow=7)

for(i in 1:5){
    bam_list_i <- bam_list[sample(11, 7, replace=F),]
    write.table(bam_list_i, row.names=F, col.names=F, quote=F)
}
