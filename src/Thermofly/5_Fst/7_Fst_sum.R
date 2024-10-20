# This code with calculate the mean, median and range for the unweighted and weighted Fst values


argv <- commandArgs(T)
file1<-argv[1]
file2<-argv[2]

fst_unweighted<-read.table(paste0(file1))
fst_weighted<-read.table(paste0(file2))

fst_unweighted_mean<-mean(fst_unweighted[,1])
fst_unweighted_median<-median(fst_unweighted[,1])
fst_unweighted_sd<-sd(fst_unweighted[,1])
fst_unweighted_range<-range(fst_unweighted[,1])

fst_unweighted_summary<-data.frame(fst_unweighted_mean, fst_unweighted_median, fst_unweighted_sd, fst_unweighted_range)
write.table(fst_unweighted_summary, "fst_unweighted_summary.txt", quote=F)

fst_weighted_mean<-mean(fst_weighted[,1])
fst_weighted_median<-median(fst_weighted[,1])
fst_weighted_sd<-sd(fst_weighted[,1])
fst_weighted_range<-range(fst_weighted[,1])

fst_weighted_summary<-data.frame(fst_weighted_mean, fst_weighted_median, fst_weighted_sd, fst_weighted_range)
write.table(fst_weighted_summary, "fst_weighted_summary.txt", quote=F)
