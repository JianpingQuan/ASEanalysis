###PCA plot of WGS individual
setwd("C:\\Users\\qjp\\Desktop\\Manuscript\\WGS\\data")

##pca value by gcta 
pca <- read.table("snp.gcta.eigenvec.txt", header = F, row.names=1)

#data inpute
library(ggplot2)
library(ggord)
pca <- pca[,-1]
colnames(pca) <- paste0("pca",1:10)
pca1_2 <- pca[,1:2]

design <- data.frame(sample = rownames(pca), group = c(rep("Duroc", 7), rep("Lulai", 8)))

Groups <- factor(design$group)

commun_theme <- theme_bw() + theme(panel.grid.major.y = element_line(linetype = "dashed",linewidth = 1), 
    panel.grid.major.x = element_line(size = 1), axis.title.x = element_text(face = "bold", size = 10), 
    axis.title.y = element_text(face = "bold", size = 10),panel.border = element_rect(linewidth = 2))

#figure output
pdf("pca_wgs.gcta.pdf", width = 5, height = 5)

ggplot(pca1_2, aes(x = pca1, y = pca2, color = Groups)) + geom_point(alpha = .7, size = 1.5) + 
    scale_color_manual(values = c("#5CB85CFF","#9632B8FF")) + labs(x = paste("PCA1 (6.1%)"), y = paste("PCA2 (1.7%)"), title = "") + stat_ellipse(level = 0.95, show.legend = F, linetype = "dashed", linewidth = 0.6) + commun_theme

dev.off()
