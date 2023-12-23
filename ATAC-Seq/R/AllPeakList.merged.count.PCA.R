setwd("F:\\Manuscript\\ATAC\\data")

##PCA Analysis
library(ggplot2)
library(dplyr)
datTPM_filter <- read.table("AllPeakList.merged.TPM2.txt", header=T, row.names = 1)

group <- read.table("Meta.txt",header = T, row.names = 1)

group <- filter(group, !(rownames(group) %in% c("LD70_M_Pl3")))
tpm <- t(datTPM_filter)
tpm1 <- tpm[,which(colSums(tpm)!=0)]  
tpm1 <- as.data.frame(tpm1)
tpm_noPl <- filter(tpm1, !grepl("Pl",rownames(tpm1)))
group_noPl <- filter(group, !grepl("Pl",rownames(group)))
stage <- c(40,70,115,168)


###All Stgae sperated and different tissues
for(i in 1:4){
  assign(paste0("tpm_",stage[i]), filter(tpm1, grepl(stage[i],rownames(tpm1))))
  assign(paste0("group_",stage[i]), filter(group, grepl(stage[i], rownames(group))))
  if(stage[i]=="40"){
    color <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF")
    shape<-c(21,22,3,24)
  } else if(stage[i]=="70"){
    color <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF")
    shape<-c(21,22,3,24)
  } else {
    color <- c("#D43F3AFF","#EEA236FF","#5CB85CFF")
    shape<-c(21,22,3)
  }
  pca = prcomp(get(paste0("tpm_",stage[i])),scale = T)
  pca_data = as.data.frame(pca$x)
  print(paste0("tpm_",stage[i]))
  print(summary(pca)$importance["Proportion of Variance", 1:2])
  pca_data = as.data.frame(cbind(pca_data,get(paste0("group_",stage[i]))))
  assign(paste0("p",stage[i]), ggplot(data = pca_data, aes(PC1,PC2)) + geom_point(size=2, aes(fill=Tissue,color=Tissue,shape=Tissue))+ stat_ellipse(aes(color = Tissue),lwd=1)+
           geom_hline(yintercept=0,linetype=4,color="grey") + 
           geom_vline(xintercept=0,linetype=4,color="grey") + 
           scale_shape_manual(values=shape) +#
           scale_fill_manual(values=color) +#
           scale_color_manual(values=color) +
           theme_classic())
  
}

pdf("PCA_40to168.pdf", height=4, width=6)
p40
p70
p115
p168
dev.off()



###TPM_Br_Li_Mu
shape<-c(21,22,3)
color <- c("#D43F3AFF","#EEA236FF","#5CB85CFF")
pca = prcomp(tpm_noPl,scale = T)
pca_data = as.data.frame(pca$x)
pca_data = as.data.frame(cbind(pca_data,group_noPl))
ggplot(data = pca_data, aes(PC1,PC2)) + geom_point(size=2, aes(fill=Tissue,color=Tissue,shape=Tissue))+ stat_ellipse(aes(color = Tissue),lwd=1)+
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  scale_shape_manual(values=shape) +#
  scale_fill_manual(values=color) +#
  scale_color_manual(values=color) +
  theme_classic()

###TPM_AllTissues
pca = prcomp(tpm1,scale = T)
pca_data = as.data.frame(pca$x)
pca_data = as.data.frame(cbind(pca_data,group))
shape<-c(21,22,3,24)
color <- c("#D43F3AFF","#EEA236FF","#5CB85CFF","#46B8DAFF")
ggplot(data = pca_data, aes(PC1,PC2)) + geom_point(size=2, aes(fill=Tissue,color=Tissue,shape=Tissue))+ stat_ellipse(aes(color = Tissue),lwd=1)+
  geom_hline(yintercept=0,linetype=4,color="grey") + 
  geom_vline(xintercept=0,linetype=4,color="grey") + 
  scale_shape_manual(values=shape) +#
  scale_fill_manual(values=color) +#
  scale_color_manual(values=color) +
  theme_classic()
  
  
