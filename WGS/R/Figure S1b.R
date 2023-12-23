###relatedness matrix heatmap and ibs matrix pca plot

#Relatedness_matrix_heatmap
setwd("C:\\Users\\qjp\\Desktop\\Manuscript\\WGS\\data")
library(pheatmap)
library(magrittr)
library(ggplot2)
cc1 = grDevices::colorRampPalette(c("darkblue","white","red") %>% ggplot2::alpha(0.7))
data <- read.table("allChr.variant.relatedness.txt", header = F)
sample <- c("DB1","DB2","DB3","DB4","DS2","DS3","DS4","LB1","LB2","LB3","LB4","LS1","LS2","LS3","LS4")
colnames(data) <- sample
rownames(data) <- colnames(data)
scale(t(data)) %>%t%>%dist%>%hclust%>%cutree(2)%>%unname -> row_hclust_cutree
row_hclust_cutree[which(row_hclust_cutree==1)]<- "Duroc"
row_hclust_cutree[which(row_hclust_cutree==2)]<- "Lulai"
annotation_row = data.frame(Breed = row_hclust_cutree)
rownames(annotation_row) <-  rownames(data)
annotation_colors = list(Breed=c(Duroc="#008B45CC",Lulai="#631879CC"))

p1 <- pheatmap::pheatmap(data,scale = "none",
                        cutree_rows=2,
                        cutree_cols = 2,
                        color=cc1(50),show_rownames = T,
                        show_colnames = F,
                        annotation_row = annotation_row, 
                        annotation_colors = annotation_colors)
ggsave("F0_Relatedness.heatmap.pdf",p1, width = 5, height = 4)
