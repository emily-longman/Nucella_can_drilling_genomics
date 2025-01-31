# Perform a pca on the covariance matrix

install.packages(c('data.table', 'ggplot2', 'ggpubr'))
library(data.table)
library(ggplot2)
library(ggpubr)

# Prior to running script, move metadata (in the format of the bam list used) and the cov matrix to script dir

# Load metadata
# Note: Make sure metadata order matches the bamlist used to produce the covariance matrix!
meta_data<-fread("Thermofly_D.basisetae.tsv", header=T)

# Load cov matrix
cov_mat <- as.matrix(read.table("Thermofly_SNPs_reduced_minInd_17_depth_6_minMaf_0.1.cov")) 
pca<-eigen(cov_mat)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

# Add column names (i.e., PCs)
nPC<-dim(pca$vectors)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca.mat)<-c(col_PC)

# Add rownames
rownames(pca.mat)<-meta_data$sampleId

# Two ways to calculate variance
# Calculate varsum(eigen_mats$values[eigen_mats$values>=0]
var1<-round(pca$values[1]*100/sum(pca$values[pca$values>=0]),3)
var2<-round(pca$values[2]*100/sum(pca$values[pca$values>=0]),3)
var3<-round(pca$values[3]*100/sum(pca$values[pca$values>=0]),3)
var4<-round(pca$values[4]*100/sum(pca$values[pca$values>=0]),3)

## How much variance is explained by the first few PCs?
var <- round(pca$values*100/sum(pca$values[pca$values>=0]),3)
var[1:5] # First 5 PCs
# Plot the eigenvalues of the PCA:
barplot(var, xlab="Eigenvalues of the PCA", ylab="Proportion of variance explained")

# Make kmeans for 2 groups on PC1
kmeans_res<-kmeans(as.matrix(pca.mat[,1]), c(min(pca.mat[,1]), median(pca.mat[,1]), max(pca.mat[,1])))
k_ss<-round(kmeans_res$betweenss/kmeans_res$totss,2)

# Combine metadata and PCs
data <- cbind(meta_data, pca.mat[,1:4])
write.table(data, "Thermofly_basisetae_PCs.csv", col.names = T, row.names = F, quote = F, sep = "\t")


#plot pca

cols=c("blue", "red")

ggscatter(data, x = "PC1", y = "PC2",
          color = "city") +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1],"%)"), y = paste0("PC2: (",var[2],"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))
ggsave("Basisetae_PC1.PC2.jpeg", width = 8, height = 6, device='jpeg', dpi=300)


ggscatter(data, x = "PC3", y = "PC4",
          color = "city") +
  theme_bw(base_size = 18, base_family = "Arial") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC3: (",var[3],"%)"), y = paste0("PC4: (",var[4],"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))
ggsave("Basisetae_PC3.PC4.jpeg", width = 8, height = 6, device='jpeg', dpi=300)
