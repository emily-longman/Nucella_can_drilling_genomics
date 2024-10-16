# Perform a pca on the covariance matrix
argv <- commandArgs(T)
INPUT <- argv[1]
BAM <- argv[2]
META <- argv[3]

install.packages(c('data.table', 'ggplot2'))
library(data.table)
library(ggplot2)

# Load metadata
meta_data<-fread(META, header=T)
meta_data_bas<-meta_data[which(meta_data$species_initial=="bas"),]

# Load cov matrix
print(paste("read cov matrix", INPUT))
cov_mat<-as.matrix(read.table(INPUT), header=F)
pca<-eigen(cov_mat)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
nPC<-dim(pca$vectors)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca.mat)<-c(col_PC)

#add rownames
rownames(pca.mat)<-meta_data_bas$sampleId
#bam_names<-read.table(BAM,header=F)
#rownames(pca.mat)<-bam_names$V1

#calculate varsum(eigen_mats$values[eigen_mats$values>=0]
var1<-round(pca$values[1]*100/sum(pca$values[pca$values>=0]),2)
var2<-round(pca$values[2]*100/sum(pca$values[pca$values>=0]),2)
var3<-round(pca$values[3]*100/sum(pca$values[pca$values>=0]),2)
var4<-round(pca$values[4]*100/sum(pca$values[pca$values>=0]),2)

#make kmeans for 2 groups on PC1
kmeans_res<-kmeans(as.matrix(pca.mat[,1]), c(min(pca.mat[,1]), median(pca.mat[,1]), max(pca.mat[,1])))
k_ss<-round(kmeans_res$betweenss/kmeans_res$totss,2)

#save 4PCS eigenvalues and k means SS
write.table(pca.mat[,1:4], paste0(INPUT,".pca"), quote=F)
write.table(c(var1,var2,var3,var4,k_ss), paste0(INPUT,".eig"), quote=F)

#plot pca

cols=c("#94D2BD", "#6d597a")

ggplot(pca.mat, x = pca.mat[,1], y = pca.mat[,2], color = "city")+
  theme_bw(base_size = 13, base_family = "Times") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1]*100,"%)"), y = paste0("PC2: (",var[2]*100,"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))

plot(pca.mat[,1], pca.mat[,2], pch=20, ylab=paste("PC2", var2), xlab=paste("PC1", var1),col=meta_data_bas$city, main=paste("k_SS",k_ss))
pca_1.2 <- plot(pca.mat[,1], pca.mat[,2], pch=20, ylab=paste("PC2", var2), xlab=paste("PC1", var1),col=kmeans_res$cluster, main=paste("k_SS",k_ss))
ggsave("pca_1.2.pdf", width = 5, height = 5, dpi=200)
#dev.off()

pca_3.4 <- plot(pca.mat[,3], pca.mat[,4], pch=20, ylab=paste("PC4", var4), xlab=paste("PC3",var3))
ggsave("pca_3.4.pdf", width = 5, height = 5, dpi=200)
#dev.off()