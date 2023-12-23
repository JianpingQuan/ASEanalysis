
##the result of tissue-specific gene expression and chromatin accessibility intensity are visualized with heat maps to see if the same cluster can be found

##
setwd("F:\\Manuscript\\ATAC\\figure\\ori_RNA_ATAC_heatmap")
library(dplyr)

GeneOpen <- read.csv("TisSpeGeneGroupbyPeakintensity_tisMean_promoter.CSV", header = T, row.names = 1) %>% arrange(.,rownames(.))

dataExp <- read.csv("F:\\ATAC\\processFiles\\RNA-tissue_specific0711\\GeneExpressionAllSamples.fittedValue.csv",header=T,row.names=1)

GeneExpFilter <- dataExp %>% filter(.,rownames(.) %in% rownames(GeneOpen)) %>% arrange(.,rownames(.)) %>% select(.,!contains(c("DL70_F_Li3","LD168_F_Br1","LD168_F_Li1","LD168_F_Mu1","LD168_F_Br2","LD168_F_Li2","LD168_F_Mu2","DL168_F_Br1","DL168_F_Li1","DL168_F_Mu1")))


##Gene expression was determined by tissue
Br <- GeneExpFilter %>% select(., contains("Br")) %>% rowMeans()
Li <- GeneExpFilter %>% select(., contains("Li")) %>% rowMeans()
Mu <- GeneExpFilter %>% select(., contains("Mu")) %>% rowMeans()
Pl <- GeneExpFilter %>% select(., contains("Pl")) %>% rowMeans()

TisCombExp <- data.frame(Brain=Br, Liver=Li, Muscle=Mu, Placenta=Pl)
rownames(TisCombExp) <- rownames(GeneExpFilter)
TisCombExp$TisInfo <- GeneOpen$TisInfo
write.csv(TisCombExp,"TisSpe_geneExp_tisMean_promoter.csv", quote = F, row.names =T)


##Plot
##RNA-seq
library(pheatmap)
library(dplyr)
setwd("F:\\Manuscript\\ATAC\\figure\\ori_RNA_ATAC_heatmap")
RNA_plot <- read.csv("TisSpe_geneExp_tisMean_promoter.csv", header =T ,row.names = 1) 

RNA_plot <- arrange(RNA_plot, TisInfo)[,1:4] %>% filter(.,!rownames(.) %in% c("ENSSSCG00000023320"))

cc1 = grDevices::colorRampPalette(c("#a6dba0","#f7f7f7","#54278f") %>% ggplot2::alpha(0.8))

annotation_row = data.frame(Type= factor(rep(c("BrSpe","LiSpe","MuSpe", "PlSpe"), c(343, 172, 125, 232))))

rownames(annotation_row) <- rownames(RNA_plot)

ann_colors = list(Type=c(BrSpe="#ca0020", LiSpe="#f4a582", MuSpe="#92c5de", PlSpe="#0571b0"))

pdf("Tissue_specific_Exp_heatmap_promoter_noScale.pdf", width=4, height=6)
exp <- pheatmap(log2(RNA_plot+0.0001),
				#scale = "row",
                cluster_rows = F,
                cluster_cols=F,
                na_col="gray",
                color = cc1(50),
                show_rownames = F,
                show_colnames = T,
                fontsize_col=10,
                angle_col=0,
                annotation_row = annotation_row,
                annotation_colors = ann_colors)	
dev.off()				


##ATAC-Seq
library(pheatmap)
library(dplyr)
setwd("F:\\Manuscript\\ATAC\\figure\\ori_RNA_ATAC_heatmap")
ATAC_plot <- read.csv("TisSpeGeneGroupbyPeakintensity_tisMean_promoter.CSV", header =T ,row.names = 1) 

ATAC_plot <- ATAC_plot %>% filter(.,!rownames(.) %in% c("ENSSSCG00000023320")) 

ATAC_plot <- arrange(ATAC_plot, TisInfo)[,1:4]

cc2 = grDevices::colorRampPalette(c("#ffffcc","#41b6c4","#253494") %>% ggplot2::alpha(0.8))

annotation_row = data.frame(Type= factor(rep(c("BrSpe","LiSpe","MuSpe", "PlSpe"), c(343, 172, 125, 232))))

rownames(annotation_row) <- rownames(RNA_plot)

ann_colors = list(Type=c(BrSpe="#ca0020", LiSpe="#f4a582", MuSpe="#92c5de", PlSpe="#0571b0"))

pdf("Tissue_specific_ATAC_heatmap_promoter_noScale.pdf", width=4, height=6)
exp <- pheatmap(log2(ATAC_plot+0.0001),
				#scale = "row",
                cluster_rows = F,
                cluster_cols=F,
                na_col="gray",
                color = cc2(50),
                show_rownames = F,
                show_colnames = T,
                fontsize_col=10,
                angle_col=0,
                annotation_row = annotation_row,
                annotation_colors = ann_colors)	
dev.off()				