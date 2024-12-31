## Extract per drilling maf in a single file both drilled and non drilled
# Take the outputs from running ANGSD on each drilling group individually and convert the .maf output into mafs for all pops as inputs for LFMM, RDA, BayPass

# This script will extract per maf in a single file for each drilling group

# ================================================================================== #

# Load packages
library(dplyr, lib="/gpfs1/home/e/l/elongman/R/x86_64-pc-linux-gnu-library/4.4")

# ================================================================================== #

# Read inputs 

argv <- commandArgs(T)
MIN_MAF <- argv[1]
PERC_IND <- argv[2]
#ANGSD_PATH<- argv[3]
MIN_DEPTH<-argv[3]
MAX_DEPTH_FACTOR<-argv[4]

# ================================================================================== #

# Read sites files for LD-pruned and filtered SNPs
sites<-read.table("sites_info/sites_all_maf_pruned", header=F)
colnames(sites)<-c("chromo", "position", "major", "minor")

# Read pop file (list of population names - i.e., Drilled, Not.Drilled)
pop<-read.table("guide_files/Phenotype_Drilling_groups.txt", header=F)
npop<-dim(pop)[1]
pop_group<-"group" #unlist(strsplit(unlist(strsplit(POP,"/"))[2],".txt"))

# Join by chromosome and position the sites and the frequencies in each population
print("join by chromosome and position the sites and the frequencies in each group/pop")
MAFall<-sites
for (i in 1:npop)
  {
    pi<-pop[i,1]
    MAFi<-read.delim(paste0("genotype_likelihoods_drilling_pruned/",pi,"_SNPs_maf",MIN_MAF,"_pctind",PERC_IND,"_mindepth",MIN_DEPTH,"_maxdepth",MAX_DEPTH_FACTOR,"_pval1e6_pruned.mafs"), header=T)
    MAFi<-MAFi[,c(1,2,6,7)]
    colnames(MAFi)<-c("chromo", "position", paste("freq", pi, sep=""),paste("n", pi, sep=""))
    head(MAFi)
    MAFall<-left_join(MAFall, MAFi, by=c("chromo", "position"))
}
head(MAFall)

#use dim to get number of rows (for snps), and non meta data columns (for pops)
nSNP<-dim(MAFall)[1]
print(paste("total nb of snp for which we ran the analysis = ", dim(MAFall)[1]))
print(paste("total nb of pop for which we ran the analysis = ", (dim(MAFall)[2]-4)/2))

which (MAFall=="NA")
#select the position which are not NA
MAFall<-MAFall[which((rowSums(MAFall[,5:dim(MAFall)[2]])>=0)),]
nSNP_no_na<-dim(MAFall)[1]
print(paste("total number of snp kept because they were covered in all populations by the chosen % of ind = ", dim(MAFall)[1]))

# ================================================================================== #

# Write the mafs in all pops/groups
write.table(MAFall, paste0("genotype_likelihoods_drilling_pruned/by_",pop_group,"_",MIN_MAF,"_pctind",PERC_IND,"_mindepth",MIN_DEPTH,"_maxdepth",MAX_DEPTH_FACTOR,".pruned.mafs"), quote=F, sep=" ")
# Write the list of SNP infered in all populations
write.table(MAFall[,1:4], paste0("genotype_likelihoods_drilling_pruned/by_",pop_group,"_",MIN_MAF,"_pctind",PERC_IND,"_mindepth",MIN_DEPTH,"_maxdepth",MAX_DEPTH_FACTOR,".pruned.snps"), quote=F, sep=" ")

# ================================================================================== #

# Write table for Baypass
#row are snp
#col are pop (2 col per pop separated by a space, with nAllele Maj, nAllele min)
print("format for baypass")
BayPass_matrix<-matrix(ncol=(npop*2), nrow=nSNP_no_na)
for (i in 1:npop)
	{
	BayPass_matrix[,(2*i-1)]<- round(MAFall[,4+(2*i-1)]*2*MAFall[,5+(2*i-1)],0)
	BayPass_matrix[,(2*i)]<- round((1-MAFall[,4+(2*i-1)])*2*MAFall[,5+(2*i-1)],0)
	}
head(BayPass_matrix)
write.table(BayPass_matrix, paste0("Outliers/baypass/by_",pop_group,"_",MIN_MAF,"_pctind",PERC_IND,"_mindepth",MIN_DEPTH,"_maxdepth",MAX_DEPTH_FACTOR,".mafs.pruned.baypass"), quote=F, row.names = F, col.names = F)

# ================================================================================== #

# Format for RDA & lfmm 
#row are snp
#col are pop
print("format for rda/lfmm")
freq_col<-seq(5,(3+npop*2), by=2)
maf_matrix<-(MAFall[,freq_col])
colnames(maf_matrix)<-pop[,1]
rownames(maf_matrix)<-paste(MAFall$chromo, MAFall$position, sep="_")
head(maf_matrix)

# Write table for lfmm
write.table(t(maf_matrix), paste0("Outliers/lfmm/by_",pop_group,"_",MIN_MAF,"_pctind",PERC_IND,".lfmm"), quote=F, row.names = F, col.names = F)


# Write table for rda add rownames
write.table(maf_matrix, paste0("Outliers/rda/by_",pop_group,"_",MIN_MAF,"_pctind",PERC_IND,".mafs.rda"), quote=F)
