
##Plot the ibs distance mtarix between RNA-seq sample according to variants calling results
setwd("C:\\Users\\qjp\\Desktop\\Manuscript\\RNA\\Ori\\data\\IBS_heatmap")

dat_now <- read.table("Autosome.allSample.ibs.matrix.txt", header = F)

meta <- read.table("meta_data3.txt", header =T, row.names =1)

rownames(dat_before) <- rownames(dat_now) <- rownames(meta)

colnames(dat_before) <- colnames(dat_now) <- rownames(meta)

library(pheatmap)

color = colorRampPalette(c("#f03b20","white"))(50)

annotation_col = data.frame(Cross = meta$Cross, Stage = meta$Stage)
rownames(annotation_col) = rownames(meta)

annoColor <- list(Cross=c(Duroc_X_Lulai="navy", Lulai_X_Duroc="firebrick3"), Stage=c(F40="#b3cde3", F70="#8c96c6", D1="#8856a7", D168="#810f7c"))

# Display row and color annotations
p <- pheatmap::pheatmap(dat_now, color=color,legend = FALSE, annotation_col = annotation_col,annotation_colors = annoColor, show_colnames = TRUE, show_rownames = FALSE) 

pdf("IBS_heatmap.pdf",width=20, height=20)
p
dev.off()
