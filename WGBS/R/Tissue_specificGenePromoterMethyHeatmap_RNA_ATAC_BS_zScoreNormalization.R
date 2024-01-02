library(dplyr)
library(pheatmap)
library(dtwclust)

Tis_spe  <- read.csv("C:\\Users\\qjp\\Desktop\\RNA-Tissue_specific.genes.csv", header = T)

##RNA-seq data
Exp <- read.csv("F:\\RNA\\RawTables\\OriginalRef\\GeneCountMatrix.tpmUnit.csv", header = T, row.names = 1) %>% select(., contains("70"))
Exp_Mean <- data.frame(BrE=rowMeans(select(Exp, contains("Br"))), LiE=rowMeans(select(Exp, contains("Li"))), MuE=rowMeans(select(Exp, contains("Mu"))), PlE=rowMeans(select(Exp, contains("Pl")))) 

#Exp_Mean_tisSpe <- filter(Exp_Mean, rownames(Exp_Mean) %in% Tis_spe$Gene)

#Exp_Mean_tisSpe$Max <- apply(Exp_Mean_tisSpe,1,function(x) {x[which.max(x)]}) %>% as.numeric()
#Exp_Mean_tisSpe_filter <- Exp_Mean_tisSpe %>% filter(.,Max <=100)
#Exp_Mean_tisSpe_filter_final <- Exp_Mean_tisSpe_filter[,1:4]

##ATAC-seq data
ATAC <- read.csv("C:\\Users\\qjp\\Desktop\\ATAC\\ATACDataUseForAnnotation.txt.GroupByGene0711.csv", header = T, row.names = 1)
ATAC_Mean <- data.frame(BrA=rowMeans(select(ATAC, contains("Br"))), LiA=rowMeans(select(ATAC, contains("Li"))), MuA=rowMeans(select(ATAC, contains("Mu"))), PlA=rowMeans(select(ATAC, contains("Pl")))) 

##BS data
sampleName <- read.table("C:\\Users\\qjp\\Desktop\\genomeFeature20220823\\data\\sampleName.txt", header = F)
Promoter_nCG <-read.table("C:\\Users\\qjp\\Desktop\\genomeFeature20220825\\data\\AllSample.gene1stExon.methy.nCG.groupby.txt")
Promoter_nX <- read.table("C:\\Users\\qjp\\Desktop\\genomeFeature20220825\\data\\AllSample.gene1stExon.methy.nX.groupby.txt")
### data preporcess
rownames(Promoter_nCG) <- Promoter_nCG[,4]
Promoter_nCG <- Promoter_nCG[,5:51]
colnames(Promoter_nCG) <- sampleName$V1
rownames(Promoter_nX) <- Promoter_nX[,4] 
Promoter_nX <- Promoter_nX[,5:51]
colnames(Promoter_nX) <- sampleName$V1
### 
MethLevel <- Promoter_nX/Promoter_nCG
MethLevel2 <- MethLevel[complete.cases(MethLevel[,1:47]),]
MethLevel_Z <- zscore(MethLevel2) %>% as.data.frame()


### Methylation level (beta value) calculate
Promoter_MethLevel_mean <- data.frame(BrM=rowMeans(select(MethLevel_Z, contains("_B"))), LiM=rowMeans(select(MethLevel_Z, contains("_L"))), MuM=rowMeans(select(MethLevel_Z, contains("_M"))), PlM=rowMeans(select(MethLevel_Z, contains("_P"))))

##
#geneCom <- intersect(rownames(Exp_Mean_tisSpe_filter_final), rownames(ATAC_Mean)) %>% intersect(., rownames(Promoter_MethLevel_mean))
geneCom <- intersect(Tis_spe$Gene,rownames(Exp_Mean)) %>% intersect(., rownames(ATAC_Mean)) %>% intersect(., rownames(Promoter_MethLevel_mean))

Exp_Mean_com <- filter(Exp_Mean, rownames(Exp_Mean) %in% geneCom) %>% arrange(.,rownames(.))
ATAC_Mean_com <- filter(ATAC_Mean, rownames(ATAC_Mean) %in% geneCom) %>% arrange(.,rownames(.))
BS_Mean_com <- filter(Promoter_MethLevel_mean, rownames(Promoter_MethLevel_mean) %in% geneCom) %>% arrange(.,rownames(.))

##
Tis_spe_com <- filter(Tis_spe, Gene %in% geneCom) %>% arrange(.,Gene)

Exp_Mean_com$Tissue <- Tis_spe_com$Tissue
ATAC_Mean_com$Tissue <- Tis_spe_com$Tissue
BS_Mean_com$Tissue <- Tis_spe_com$Tissue

Exp_Mean_com_sort <- Exp_Mean_com %>% arrange(., Tissue)
ATAC_Mean_com_sort <- ATAC_Mean_com %>% arrange(., Tissue)
BS_Mean_com_sort <- BS_Mean_com %>% arrange(., Tissue)

Exp_plot <- log2(Exp_Mean_com_sort[,1:4]+0.0001)
ATAC_plot <- log2(ATAC_Mean_com_sort[,1:4]+0.0001)
BS_plot <- BS_Mean_com_sort[,1:4]

cc1 = grDevices::colorRampPalette(c("#f03b20","#fecc5c","#ffffb2") %>% ggplot2::alpha(0.8))
cc2 = grDevices::colorRampPalette(c("#ffffcc","#41b6c4","#253494") %>% ggplot2::alpha(0.8))
cc3 = grDevices::colorRampPalette(c("#f2f0f7","#9e9ac8","#54278f") %>% ggplot2::alpha(0.8))


#annotation_row = data.frame(Type= factor(rep(c("BrSpe","LiSpe","MuSpe", "PlSpe"), c(592, 287, 199, 231))))
##promoter 
#annotation_row = data.frame(Type= factor(rep(c("BrSpe","LiSpe","MuSpe", "PlSpe"), c(321, 154, 108, 201))))
annotation_row = data.frame(Type= factor(rep(c("BrSpe","LiSpe","MuSpe", "PlSpe"), 
                                             c(204, 114, 76, 132))))
rownames(annotation_row) <- rownames(Exp_plot)


ann_colors = list(Type=c(BrSpe="#ca0020", LiSpe="#f4a582", MuSpe="#92c5de", PlSpe="#0571b0"))


exp <- pheatmap(Exp_plot,
                cluster_rows = F,
                cluster_cols=F,
                na_col="gray",
                color = cc3(50),
                show_rownames = F,
                show_colnames = T,
                fontsize_col=8,
                angle_col=90,
                annotation_row = annotation_row,
                annotation_colors = ann_colors
)


atac <- pheatmap(ATAC_plot,
                 cluster_rows = F,
                 cluster_cols=F,
                 na_col="gray",
                 color = cc2(50),
                 show_rownames = F,
                 show_colnames = T,
                 fontsize_col=8,
                 angle_col=90,
                 annotation_row = annotation_row,
                 annotation_colors = ann_colors
)

bs <- pheatmap(BS_plot,
               cluster_rows = F,
               cluster_cols=F,
               #gaps_col = c(1,2,3),
               #gaps_row = c(363,537,663),
               na_col="gray",
               color = cc1(50),
               show_rownames = F,
               show_colnames = T,
               fontsize_col=8,
               angle_col=90,
               annotation_row = annotation_row,
               annotation_colors = ann_colors
)
