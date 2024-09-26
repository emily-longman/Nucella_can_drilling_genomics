argv <- commandArgs(T)
INFILE <- argv[1]
OUTFILE_sites <- argv[2]
OUTFILE_regions <- argv[3]

maf<-read.table(INFILE, header=T)
head (maf)
sites<-maf[,1:4]
sites_order <- sites[order(sites$chromo),] 
sites_order$chromo <-as.factor(sites_order$chromo)
regions_order<-levels(sites_order$chromo)
write.table(sites_order, OUTFILE_sites, row.names=F, col.names=F, sep="\t", quote=F)
write.table(regions_order, OUTFILE_regions, row.names=F, col.names=F, sep="\t", quote=F)