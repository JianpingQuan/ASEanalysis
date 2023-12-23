setwd("F:\\Manuscript\\ATAC\\data")
library(dplyr)
peak_gene <- read.csv("AllPeakList.merged.ReadCount.PeakToGene_Feature.csv", header = T)

TisSpeGene <- read.csv("F:\\ATAC\\processFiles\\RNA-tissue_specific0711\\RNA-Tissue_specific.genes.csv", header = T) %>% arrange(.,Var1)

##extract the peaks that can be annotated to tis-spe gene
TisSpePeak_gene <- filter(peak_gene, gene %in% TisSpeGene$Var1) %>% arrange(.,Peaks)

##TPM
PeakTPM <- read.table("F:\\Manuscript\\ATAC\\data\\AllPeakList.merged.TPM.txt", header = T, row.names=1)

##Read COunt
PeakRC <- read.table("F:\\Manuscript\\ATAC\\data\\AllPeakList.merged.ReadCount.txt", header = T, row.names=1)

##TPM filtering 
PeakTPMFilter <- filter(PeakTPM, rownames(PeakTPM) %in% TisSpePeak_gene$Peaks) %>% select(.,!contains(c("DL70_F_Li3","LD168_F_Br1","LD168_F_Li1","LD168_F_Mu1","LD168_F_Br2","LD168_F_Li2","LD168_F_Mu2","DL168_F_Br1","DL168_F_Li1","DL168_F_Mu1"))) %>% arrange(.,rownames(.))

##ReadCount filtering
PeakRCFilter <- filter(PeakRC, rownames(PeakRC) %in% TisSpePeak_gene$Peaks) %>% select(.,!contains(c("DL70_F_Li3","LD168_F_Br1","LD168_F_Li1","LD168_F_Mu1","LD168_F_Br2","LD168_F_Li2","LD168_F_Mu2","DL168_F_Br1","DL168_F_Li1","DL168_F_Mu1"))) %>% arrange(.,rownames(.))

##add gene name
PeakTPMFilter$Gene  <- TisSpePeak_gene$gene
PeakRCFilter$Gene <- TisSpePeak_gene$gene


##merge the intensity of peaks from same gene
PeakTPMFilterGroupBy <- aggregate(. ~ Gene, data=PeakTPMFilter, FUN=sum)

##
write.csv(PeakTPMFilterGroupBy, "TisSpeGeneGroupbyPeakintensity.CSV", quote = F, row.names = FALSE)


##Find the average value of each tissue
Br <- PeakTPMFilterGroupBy %>% select(., contains("Br")) %>% rowMeans()
Li <- PeakTPMFilterGroupBy %>% select(., contains("Li")) %>% rowMeans()
Mu <- PeakTPMFilterGroupBy %>% select(., contains("Mu")) %>% rowMeans()
Pl <- PeakTPMFilterGroupBy %>% select(., contains("Pl")) %>% rowMeans()

TisComb <- data.frame(Brain=Br, Liver=Li, Muscle=Mu, Placenta=Pl)
rownames(TisComb) <- PeakTPMFilterGroupBy$Gene

##sort according to gene name
TisComb <- arrange(TisComb, rownames(TisComb))

TisSpeGene_com <- filter(TisSpeGene, Var1 %in% rownames(TisComb))
TisComb$TisInfo <- TisSpeGene_com$Tissue

##
write.csv(TisComb, "TisSpeGeneGroupbyPeakintensity_tisMean.CSV", quote = F, row.names =TRUE)
