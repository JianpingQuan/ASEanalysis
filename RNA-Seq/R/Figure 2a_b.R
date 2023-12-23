##1) convert reads count based gene expression to TPM based
library(dplyr)
library("DGEobj.utils")

setwd("C:\\Users\\qjp\\Desktop\\Manuscript\\RNA\\Ori\\data\\featureCount\\")

data <- read.csv("GeneExpressionAllSamples.csv", header = T, row.names = 1) %>% as.matrix()

length <- read.csv("GeneLength.csv", row.names=1)

tpm <- convertCounts(data, unit = "TPM", geneLength = length$length)

write.csv(tpm, "GeneExpressionAllSamples.tpm.csv", quote = F, row.names = TRUE)



##2) PCA analysis
library("ggplot2")
library("factoextra")
library("ggrepel")
library(dplyr)
library(FactoMineR)

setwd("C:\\Users\\qjp\\Desktop\\Manuscript\\RNA\\Ori\\data\\PCA\\featureCount_result")

tpm <- read.csv("GeneExpressionAllSamples.tpmFilterSample.csv",header = T,row.names = 1)
group <- read.table("group.txt",header = T,sep = "\t",row.names = 1)

tpm <- t(tpm) 
tpm1 <- tpm[,which(colSums(tpm)!=0)] 
dim(tpm)
dim(tpm1) 
pca <- prcomp(tpm1,scale = T) 
res.pca <- PCA(tpm1, graph=FALSE)
head(res.pca$eig)

color <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF") 
shape<-c(21,22,3,24)
legend_title = "PCA"
xlab=paste("PC1 (23.8%)") 
ylab2=paste("PC2 (9.7%)")
ylab3=paste("PC3 (7.1%)")
pca_data<-as.data.frame(pca$x)
group$age<-as.factor(group$Age)
pca_data<-as.data.frame(cbind(pca_data,group))


pdf("pca_PC1_PC2.tpm.pdf",h=5,w=5.5)
ggplot(pca_data,aes(PC1,PC2)) + 
     geom_point(size=2,aes(fill=Tissue,color=Tissue,shape=Age)) + 
     labs(x=xlab,y=ylab2) + 
     stat_ellipse(aes(color = Tissue),linewidth=1)+
     geom_hline(yintercept=0,linetype=4,color="grey") + 
     geom_vline(xintercept=0,linetype=4,color="grey") + 
     scale_shape_manual(values=shape) +
     scale_fill_manual(values=color) +
     scale_color_manual(values=color) +
     theme_classic()
 #theme(legend.position="none")
dev.off()


pdf("pca_PC1_PC3.tpm.pdf",h=5,w=5.5)
ggplot(pca_data,aes(PC1,PC3)) +
  geom_point(size=2,aes(fill=Tissue,color=Tissue,shape=Age)) + 
  labs(x=xlab,y=ylab3) + 
  stat_ellipse(aes(color = Tissue),linewidth=1)+
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  scale_shape_manual(values=shape) +
  scale_fill_manual(values=color) +
  scale_color_manual(values=color) +
  theme_classic()
#theme(legend.position="none")
dev.off()
