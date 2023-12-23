#
library(dplyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(reshape2)

setwd("F:\\Manuscript\\WGBS")

##peaks containing CpG
datBS <- read.csv("AllSamples.CpG.removeSNP.CovFilter.MethyLevel.csv", header = T, row.names = 1)

##peaks intensity (TPM)
datATAC <- read.table("F:\\Manuscript\\ATAC\\data\\AllPeakList.merged.TPM.txt", header = T, row.names=1)
datATAC <- select(datATAC, contains("70"))

##common peaks
peakCom <- intersect(rownames(datATAC), rownames(datBS))

datATAC_peakCom <- filter(datATAC, rownames(datATAC) %in% peakCom) %>% arrange(.,rownames(.))

datBS_peakCom <- filter(datBS, rownames(datBS) %in% peakCom) %>% arrange(.,rownames(.))

##
datATAC_peakComMean <- data.frame(ATACBr=rowMeans(select(datATAC_peakCom, contains("_B"))), ATACLi=rowMeans(select(datATAC_peakCom, contains("_L"))), ATACMu=rowMeans(select(datATAC_peakCom, contains("_M"))), ATACPl=rowMeans(select(datATAC_peakCom, contains("_P")))) %>% set_rownames(rownames(datATAC_peakCom)) %>% arrange(., rownames(.))

datBS_peakComMean <- data.frame(MethBr=rowMeans(select(datBS_peakCom, contains("_B"))), MethLi=rowMeans(select(datBS_peakCom, contains("_L"))), MethMu=rowMeans(select(datBS_peakCom, contains("_M"))), MethPl=rowMeans(select(datBS_peakCom, contains("_P")))) %>% set_rownames(rownames(datBS_peakCom)) %>% arrange(., rownames(.))


dataCombine <- cbind(datATAC_peakComMean, datBS_peakComMean)


##only consider the Peak that located in promoter region of genes
annotate <- read.table("F:\\Manuscript\\ATAC\\data\\PromoterPeaktoGene.txt", header=T)

dataCombine_f <- filter(dataCombine, rownames(dataCombine) %in% annotate$Peaks )

cor.test(dataCombine_f$ATACBr,dataCombine_f$MethBr, method = "spearman")
cor.test(dataCombine_f$ATACLi,dataCombine_f$MethLi, method = "spearman")
cor.test(dataCombine_f$ATACMu,dataCombine_f$MethMu, method = "spearman")
cor.test(dataCombine_f$ATACPl,dataCombine_f$MethPl, method = "spearman")



Br <- ggplot(dataCombine_f, aes(x=MethBr, y=log2(ATACBr+1)))+ 
  geom_point(alpha = 0.03, color = "#3F3F3F")+
  xlab(label = NULL)+ylab(label = "Chromatin accessibility") + 
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1.0), labels = c(0,0.25,0.5,0.75,1.0))+
  geom_smooth(color = "#EEA236", fill="gray")+theme_bw()+
  ylab(label = "Chromatin accessibility") + 
  labs(title = "Brain")+ theme(plot.title = element_text(hjust = 0.5))+
  annotate("text",x=0.75, y=10, label = "rho = -0.64", size =6) + theme(plot.title = element_text(size = 15), axis.title.y = element_text(size = 15), axis.text = element_text(size = 15))

Li <- ggplot(dataCombine_f, aes(x=MethLi, y=log2(ATACLi+1)))+ 
  geom_point(alpha = 0.03, color = "#3F3F3F")+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1.0), labels = c(0,0.25,0.5,0.75,1.0))+
  xlab(label = NULL) + ylab(label = NULL) + 
  geom_smooth(color = "#EEA236", fill="gray")+theme_bw()+
  labs(title = "Liver")+theme(plot.title = element_text(hjust = 0.5))+
  annotate("text",x=0.75, y=10, label = "rho = -0.49", size =6) + theme(plot.title = element_text(size = 15), axis.title.y = element_text(size = 15), axis.text = element_text(size = 15))

Mu <- ggplot(dataCombine_f, aes(x=MethMu, y=log2(ATACMu+1)))+ 
  geom_point(alpha = 0.03, color = "#3F3F3F")+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1.0),labels = c(0,0.25,0.5,0.75,1.0))+
  geom_smooth(color = "#EEA236", fill="gray")+theme_bw()+
  xlab(label = "Methylation level")+ylab(label = "Chromatin accessibility") + 
  labs(title = "Muscle")+theme(plot.title = element_text(hjust = 0.5))+
  annotate("text",x=0.75, y=10, label = "rho = -0.61", size =6) + theme(plot.title = element_text(size = 15), axis.title.y = element_text(size = 15),axis.title.x = element_text(size = 15), axis.text = element_text(size = 15))

Pl <- ggplot(dataCombine_f, aes(x=MethPl, y=log2(ATACPl+1)))+ 
  geom_point(alpha = 0.03, color = "#3F3F3F")+
  scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1.0), labels = c(0,0.25,0.5,0.75,1.0))+
  geom_smooth(color = "#EEA236", fill="gray")+theme_bw()+
  xlab(label = "Methylation level")+ ylab(label = NULL) + 
  labs(title = "Placenta")+theme(plot.title = element_text(hjust = 0.5))+
  annotate("text",x=0.75, y=10, label = "rho = -0.46", size =6) + theme(plot.title = element_text(size = 15), axis.title.y = element_text(size = 15),axis.title.x = element_text(size = 15), axis.text = element_text(size = 15))

all <- ggarrange(Br,Li,Mu,Pl, nrow = 2, ncol = 2)

pdf("RelationShipeBetween_ATAC_Meth_onlyPromoter.pdf", width = 7, height = 6)
all
dev.off()


png("RelationShipeBetween_ATAC_Meth.png", width = 700, height = 600, res = 600)
all
dev.off()


jpeg("RelationShipeBetween_ATAC_Meth.jpg", width = 1400, height = 1200, quality = 100)
all
dev.off()