Input="F:\\RNA_20221101\\02_ModifiedGenome\\A_GeneExpressionMatrix\\"

library(dplyr)
library(edgeR)
library(ggplot2)
library(stringr)

data_filter <- read.csv(paste0(Input,"AutosomeGeneCount_filtered2.csv"), header = T, row.names = 1)

POEgeneList <- read.csv("F:\\RNA_20221101\\02_ModifiedGenome\\C_POEgeneSelection\\result\\POE.list2.csv", header = T, row.names = 1)

data_POE <- filter(data_filter, rownames(data_filter) %in% POEgeneList$all_gene_unique)

dat_POE_Mat <- select(data_POE, contains("Mat"))
dat_POE_Pat <- select(data_POE, contains("Pat"))

Mat_addPat <- dat_POE_Mat+dat_POE_Pat+0.0002

MatBiasScore <- (dat_POE_Mat+0.0001)/Mat_addPat

colnames(MatBiasScore) <- str_sub(colnames(dat_POE_Mat),1,-7)


##The replicate bias score was averaged for the same organization over the same period
F40_Br <- select(MatBiasScore, contains("40")&contains("Br")) %>% rowMeans() %>% as.data.frame() 
F40_Li <- select(MatBiasScore, contains("40")&contains("Li")) %>% rowMeans() %>% as.data.frame()
F40_Mu <- select(MatBiasScore, contains("40")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
F40_Pl <- select(MatBiasScore, contains("40")&contains("Pl")) %>% rowMeans() %>% as.data.frame()
F70_Br <- select(MatBiasScore, contains("70")&contains("Br")) %>% rowMeans() %>% as.data.frame()
F70_Li <- select(MatBiasScore, contains("70")&contains("Li")) %>% rowMeans() %>% as.data.frame()
F70_Mu <- select(MatBiasScore, contains("70")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
F70_Pl <- select(MatBiasScore, contains("70")&contains("Pl")) %>% rowMeans() %>% as.data.frame()
D1_Br <-  select(MatBiasScore, contains("115")&contains("Br")) %>% rowMeans() %>% as.data.frame()
D1_Li <-  select(MatBiasScore, contains("115")&contains("Li")) %>% rowMeans() %>% as.data.frame()
D1_Mu <-  select(MatBiasScore, contains("115")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
D168_Br <- select(MatBiasScore, contains("168")&contains("Br")) %>% rowMeans() %>% as.data.frame()
D168_Li <- select(MatBiasScore, contains("168")&contains("Li")) %>% rowMeans() %>% as.data.frame()
D168_Mu <- select(MatBiasScore, contains("168")&contains("Mu")) %>% rowMeans() %>% as.data.frame()


##
F40_Br_im <- filter(F40_Br,.>=0.95|.<=0.05)
F40_Br_im$Group <- rep("F40_Br",times=nrow(F40_Br_im))
F40_Br_im$Gene <- rownames(F40_Br_im)

F40_Li_im <- filter(F40_Li,.>=0.95|.<=0.05)
F40_Li_im$Group <- rep("F40_Li",times=nrow(F40_Li_im))
F40_Li_im$Gene <- rownames(F40_Li_im)

F40_Mu_im <- filter(F40_Mu,.>=0.95|.<=0.05)
F40_Mu_im$Group <- rep("F40_Mu",times=nrow(F40_Mu_im))
F40_Mu_im$Gene <- rownames(F40_Mu_im)

F40_Pl_im <- filter(F40_Pl,.>=0.95|.<=0.05)
F40_Pl_im$Group <- rep("F40_Pl",times=nrow(F40_Pl_im))
F40_Pl_im$Gene <- rownames(F40_Pl_im)

F70_Br_im <- filter(F70_Br,.>=0.95|.<=0.05)
F70_Br_im$Group <- rep("F70_Br",times=nrow(F70_Br_im))
F70_Br_im$Gene <- rownames(F70_Br_im)

F70_Li_im <- filter(F70_Li,.>=0.95|.<=0.05)
F70_Li_im$Group <- rep("F70_Li",times=nrow(F70_Li_im))
F70_Li_im$Gene <- rownames(F70_Li_im)

F70_Mu_im <- filter(F70_Mu,.>=0.95|.<=0.05)
F70_Mu_im$Group <- rep("F70_Mu",times=nrow(F70_Mu_im))
F70_Mu_im$Gene <- rownames(F70_Mu_im)

F70_Pl_im <- filter(F70_Pl,.>=0.95|.<=0.05)
F70_Pl_im$Group <- rep("F70_Pl",times=nrow(F70_Pl_im))
F70_Pl_im$Gene <- rownames(F70_Pl_im)

D1_Br_im <- filter(D1_Br,.>=0.95|.<=0.05)
D1_Br_im$Group <- rep("D1_Br",times=nrow(D1_Br_im))
D1_Br_im$Gene <- rownames(D1_Br_im)

D1_Li_im <- filter(D1_Li,.>=0.95|.<=0.05)
D1_Li_im$Group <- rep("D1_Li",times=nrow(D1_Li_im))
D1_Li_im$Gene <- rownames(D1_Li_im)

D1_Mu_im <- filter(D1_Mu,.>=0.95|.<=0.05)
D1_Mu_im$Group <- rep("D1_Mu",times=nrow(D1_Mu_im))
D1_Mu_im$Gene <- rownames(D1_Mu_im)

D168_Br_im <- filter(D168_Br,.>=0.95|.<=0.05)
D168_Br_im$Group <- rep("D168_Br",times=nrow(D168_Br_im))
D168_Br_im$Gene <- rownames(D168_Br_im)

D168_Li_im <- filter(D168_Li,.>=0.95|.<=0.05)
D168_Li_im$Group <- rep("D168_Li",times=nrow(D168_Li_im))
D168_Li_im$Gene <- rownames(D168_Li_im)

D168_Mu_im <- filter(D168_Mu,.>=0.95|.<=0.05)
D168_Mu_im$Group <- rep("D168_Mu",times=nrow(D168_Mu_im))
D168_Mu_im$Gene <- rownames(D168_Mu_im)

all <- rbind(F40_Br_im,F40_Li_im,F40_Mu_im,F40_Pl_im,F70_Br_im,F70_Li_im,F70_Mu_im,F70_Pl_im,D1_Br_im,D1_Li_im,D1_Mu_im,D168_Br_im,D168_Li_im,D168_Mu_im)

imGene <- unique(sort(all$Gene))

proveIM <- read.csv("F:\\Manuscript\\正反交论文\\figures\\ImprintGenesDatabase.csv", header = T, row.names = 1)
proveIMList <- proveIM$Ensembl

##novel imprinting gene list2
noveIm <- setdiff(imGene,proveIMList)


##Bias score of each novel imprinting gene
F40_Br_novelIm <- F40_Br %>% filter(., rownames(F40_Br) %in% noveIm)
F70_Br_novelIm <- F70_Br %>% filter(., rownames(F70_Br) %in% noveIm)
D1_Br_novelIm <- D1_Br %>% filter(., rownames(D1_Br) %in% noveIm)
D168_Br_novelIm <- D168_Br %>% filter(., rownames(D168_Br) %in% noveIm)

F40_Li_novelIm <- F40_Li %>% filter(., rownames(F40_Li) %in% noveIm)
F70_Li_novelIm <- F70_Li %>% filter(., rownames(F70_Li) %in% noveIm)
D1_Li_novelIm <- D1_Li %>% filter(., rownames(D1_Li) %in% noveIm)
D168_Li_novelIm <- D168_Li %>% filter(., rownames(D168_Li) %in% noveIm)

F40_Mu_novelIm <- F40_Mu %>% filter(., rownames(F40_Mu) %in% noveIm)
F70_Mu_novelIm <- F70_Mu %>% filter(., rownames(F70_Mu) %in% noveIm)
D1_Mu_novelIm <- D1_Mu %>% filter(., rownames(D1_Mu) %in% noveIm)
D168_Mu_novelIm <- D168_Mu %>% filter(., rownames(D168_Mu) %in% noveIm)

F40_Pl_novelIm <- F40_Pl %>% filter(., rownames(F40_Pl) %in% noveIm)
F70_Pl_novelIm <- F70_Pl %>% filter(., rownames(F70_Pl) %in% noveIm)

allNovelImBS <- cbind(F40_Br_novelIm,F40_Li_novelIm,F40_Mu_novelIm,F40_Pl_novelIm,F70_Br_novelIm,F70_Li_novelIm,F70_Mu_novelIm,F70_Pl_novelIm,D1_Br_novelIm,D1_Li_novelIm,D1_Mu_novelIm,D168_Br_novelIm,D168_Li_novelIm,D168_Mu_novelIm)
colnames(allNovelImBS) <-  c("F40_Br","F40_Li","F40_Mu","F40_Pl","F70_Br","F70_Li","F70_Mu","F70_Pl","D1_Br","D1_Li","D1_Mu","D168_Br","D168_Li","D168_Mu")
write.csv(allNovelImBS, "F://Manuscript//正反交论文//figures//NovelImprintingGene_biasScore.csv", quote = F, row.names = T)


##
data_novel_POE <- filter(data_filter, rownames(data_filter) %in% noveIm)

F40_Br_novel <- select(data_novel_POE, contains("40")&contains("Br")) %>% rowMeans() %>% as.data.frame() 
F40_Li_novel <- select(data_novel_POE, contains("40")&contains("Li")) %>% rowMeans() %>% as.data.frame()
F40_Mu_novel <- select(data_novel_POE, contains("40")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
F40_Pl_novel <- select(data_novel_POE, contains("40")&contains("Pl")) %>% rowMeans() %>% as.data.frame()
F70_Br_novel <- select(data_novel_POE, contains("70")&contains("Br")) %>% rowMeans() %>% as.data.frame()
F70_Li_novel <- select(data_novel_POE, contains("70")&contains("Li")) %>% rowMeans() %>% as.data.frame()
F70_Mu_novel <- select(data_novel_POE, contains("70")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
F70_Pl_novel <- select(data_novel_POE, contains("70")&contains("Pl")) %>% rowMeans() %>% as.data.frame()
D1_Br_novel <-  select(data_novel_POE, contains("115")&contains("Br")) %>% rowMeans() %>% as.data.frame()
D1_Li_novel <-  select(data_novel_POE, contains("115")&contains("Li")) %>% rowMeans() %>% as.data.frame()
D1_Mu_novel <-  select(data_novel_POE, contains("115")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
D168_Br_novel <- select(data_novel_POE, contains("168")&contains("Br")) %>% rowMeans() %>% as.data.frame()
D168_Li_novel <- select(data_novel_POE, contains("168")&contains("Li")) %>% rowMeans() %>% as.data.frame()
D168_Mu_novel <- select(data_novel_POE, contains("168")&contains("Mu")) %>% rowMeans() %>% as.data.frame()

all_novel <- cbind(F40_Br_novel,F40_Li_novel,F40_Mu_novel,F40_Pl_novel,F70_Br_novel,F70_Li_novel,F70_Mu_novel,F70_Pl_novel,D1_Br_novel,D1_Li_novel,D1_Mu_novel,D168_Br_novel, D168_Li_novel, D168_Mu_novel)
colnames(all_novel) <- c("F40_Br","F40_Li","F40_Mu","F40_Pl","F70_Br","F70_Li","F70_Mu","F70_Pl","D1_Br","D1_Li","D1_Mu","D168_Br","D168_Li","D168_Mu")
all_novel <- log2(all_novel+0.0001)
all_novel$id <- rownames(all_novel)


all_novel_bias <- MatBiasScore %>% filter(.,rownames(.) %in% noveIm)

F40_Br <- select(all_novel_bias, contains("40")&contains("Br")) %>% rowMeans() %>% as.data.frame() 
F40_Li <- select(all_novel_bias, contains("40")&contains("Li")) %>% rowMeans() %>% as.data.frame()
F40_Mu <- select(all_novel_bias, contains("40")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
F40_Pl <- select(all_novel_bias, contains("40")&contains("Pl")) %>% rowMeans() %>% as.data.frame()
F70_Br <- select(all_novel_bias, contains("70")&contains("Br")) %>% rowMeans() %>% as.data.frame()
F70_Li <- select(all_novel_bias, contains("70")&contains("Li")) %>% rowMeans() %>% as.data.frame()
F70_Mu <- select(all_novel_bias, contains("70")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
F70_Pl <- select(all_novel_bias, contains("70")&contains("Pl")) %>% rowMeans() %>% as.data.frame()
D1_Br <-  select(all_novel_bias, contains("115")&contains("Br")) %>% rowMeans() %>% as.data.frame()
D1_Li <-  select(all_novel_bias, contains("115")&contains("Li")) %>% rowMeans() %>% as.data.frame()
D1_Mu <-  select(all_novel_bias, contains("115")&contains("Mu")) %>% rowMeans() %>% as.data.frame()
D168_Br <- select(all_novel_bias, contains("168")&contains("Br")) %>% rowMeans() %>% as.data.frame()
D168_Li <- select(all_novel_bias, contains("168")&contains("Li")) %>% rowMeans() %>% as.data.frame()
D168_Mu <- select(all_novel_bias, contains("168")&contains("Mu")) %>% rowMeans() %>% as.data.frame()

all_novel_biasMean <- cbind(F40_Br,F40_Li,F40_Mu,F40_Pl,F70_Br,F70_Li,F70_Mu,F70_Pl,D1_Br,D1_Li,D1_Mu,D168_Br, D168_Li, D168_Mu)
colnames(all_novel_biasMean) <- c("F40_Br","F40_Li","F40_Mu","F40_Pl","F70_Br","F70_Li","F70_Mu","F70_Pl","D1_Br","D1_Li","D1_Mu","D168_Br","D168_Li","D168_Mu")
all_novel_biasMean$id <- rownames(all_novel_biasMean)


library(reshape2)
all_novel_long <- melt(all_novel, value.name = c("Abundance"), id.vars = c("id"))
all_novel_biasMean_long <- melt(all_novel_biasMean, value.name = c("Bias"), id.vars = c("id"))

m3 <- cbind(all_novel_long, all_novel_biasMean_long)
m4<-m3[,c(1:3,6)]
m4 <- arrange(m4, id)

p2<-ggplot(m4,aes(variable,id,size=`Abundance`))+
geom_point(shape=21,aes(fill=Bias), position =position_dodge(0))+
theme_minimal()+xlab(NULL)+ylab(NULL)+
scale_fill_gradientn(colours=c("#00468B","#F8E5C6","#ED0000"),guide="legend")+
theme(legend.position = "bottom",legend.box = "vertical", panel.grid.major =element_blank(), legend.margin=margin(t= 0, unit='cm'), legend.spacing = unit(0,"in"))+
scale_size_continuous(range=c(-4,5)) + theme(axis.text.x=element_text(angle=45)) + theme(text = element_text(size = 18))

ggsave(file="NovelImprinting.pdf",height=8,width = 10)