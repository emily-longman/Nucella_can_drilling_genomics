# Perform a pca on the covariance matrix

argv <- commandArgs(T)
INPUT <- argv[1]
BAM <- argv[2]
META <- argv[3]

install.packages(c('data.table', 'ggplot2'))
library(data.table)
library(ggplot2)
library(ggpubr)

# Load metadata
meta_data<-fread("Thermofly_metadata.tsv", header=T)
meta_data_bas<-meta_data[which(meta_data$species_initial=="bas"),]
# Remove metadata for individuals with low coverage
meta_data_bas_reduced <- meta_data_bas[-which(meta_data$sampleId =="D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_9.1"| 
                                        meta_data$sampleId == "D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_12.1" |
                                        meta_data$sampleId == "D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_13.1" |
                                        meta_data$sampleId == "D_bas.wild.US-HI-Ola.27_1_2023.w.DP_B_O_15.1")]

# Load cov matrix
cov_mat <- as.matrix(read.table("Thermofly_SNPs_reduced_minInd_16_depth_4.cov")) 
pca<-eigen(cov_mat)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

# Add column names (i.e., PCs)
nPC<-dim(pca$vectors)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca.mat)<-c(col_PC)

# Add rownames
rownames(pca.mat)<-meta_data_bas_reduced$sampleId

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
data <- cbind(meta_data_bas_reduced, pca.mat[,1:4])
write.table(data, "Thermofly_basisetae_PCs.csv", col.names = T, row.names = F, quote = F, sep = "\t")

#plot pca

cols=c("#94D2BD", "#6d597a")
cols.all=c("pink", "red", "orange", "yellow", "green", "darkolivegreen3","darkgreen", "cyan", "skyblue",
           "blue", "darkblue", "purple", "darkorchid4", "darkkhaki", "grey", "brown", "black", "tan")

ggscatter(data, x = "PC1", y = "PC2",
          color = "city") +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1],"%)"), y = paste0("PC2: (",var[2],"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))
ggsave("Basisetae_PC1.PC2.jpeg", width = 8, height = 6, device='jpeg', dpi=300)

ggscatter(data, x = "PC1", y = "PC2",
          color = "sampleId") +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC1: (",var[1],"%)"), y = paste0("PC2: (",var[2],"%)")) +
  scale_color_manual(values=c(cols.all), name="Source population") +
  guides(colour = guide_legend(nrow = 6))

ggscatter(data, x = "PC3", y = "PC4",
          color = "city") +
  theme_bw(base_size = 13, base_family = "Arial") +
  theme(panel.background = element_blank(), 
        legend.background = element_blank(), 
        panel.grid = element_blank(), 
        plot.background = element_blank(), 
        legend.text=element_text(size=rel(.7)), 
        axis.text = element_text(size=13), 
        legend.position = "bottom") +
  labs(x = paste0("PC3: (",var[1],"%)"), y = paste0("PC4: (",var[2],"%)")) +
  scale_color_manual(values=c(cols), name="Source population") +
  guides(colour = guide_legend(nrow = 2))
ggsave("Basisetae_PC3.PC4.jpeg", width = 8, height = 6, device='jpeg', dpi=300)
