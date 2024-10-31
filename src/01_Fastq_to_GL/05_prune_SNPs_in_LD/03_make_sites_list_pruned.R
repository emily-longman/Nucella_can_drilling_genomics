# This workflow will extract the pruned SNPs to make a list of sites for subsequent analysis in angsd

argv <- commandArgs(T)
INPUT_plink <- argv[1]
INPUT_angsd <- argv[2]

lib.EKL<-/gpfs1/home/e/l/elongman/R/x86_64-pc-linux-gnu-library/4.4/00LOCK-dplyr/00new/dplyr/libs

# Install dplyr
install.packages("dplyr",lib=lib.SQL)
#Load dplyr
library("dplyr")

pruned<-read.table(INPUT_plink)
head(pruned)
colnames(pruned)<-"LG_pos"
sites<-read.table(INPUT_angsd)
head(sites)
sites$LG_pos<-paste0(sites[,1],"_", sites[,2])


sites_pruned<-inner_join(sites, pruned)
print(paste("there is a total of ", dim(sites)[1], "sites"))
print(paste("we keep ", dim(sites_pruned)[1], "pruned sites"))

write.table(sites_pruned[,1:4], paste0(INPUT_angsd,"_pruned"), col.names=F, row.names=F, quote=F, sep="\t")