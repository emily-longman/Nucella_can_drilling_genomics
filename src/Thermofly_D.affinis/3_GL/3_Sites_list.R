argv <- commandArgs(T)
infile <- argv[1]
outfile_sites <- argv[2]
outfile_regions <- argv[3]

maf<-read.table(infile, header=T)
head (maf)
sites<-maf[,1:4]
sites_order <- sites[order(sites$chromo),] 
sites_order$chromo <-as.factor(sites_order$chromo)
regions_order<-levels(sites_order$chromo)
write.table(sites_order, outfile_sites, row.names=F, col.names=F, sep="\t", quote=F)
write.table(regions_order, outfile_regions, row.names=F, col.names=F, sep="\t", quote=F)