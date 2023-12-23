library(dplyr)

OriExp <- read.csv("C:\\Users\\qjp\\Desktop\\Manuscript\\RNA\\Ori\\data\\Tissue_Stage_specific\\GeneExpressionAllSamples.fittedValue.csv", header = T, row.names = 1)

OriExp70 <- select(OriExp, contains("70"))

gene <- c("ENSSSCG00000040410","ENSSSCG00000046943","ENSSSCG00000047414","ENSSSCG00000035805","ENSSSCG00000048717")

OriExp70gene <- filter(OriExp70, rownames(OriExp70) %in% gene)

Brain <- OriExp70gene %>% select(.,contains("Br")) %>% t() %>% as.data.frame()
Liver <- OriExp70gene %>% select(.,contains("Li")) %>% t() %>% as.data.frame()
Muscle <- OriExp70gene %>% select(.,contains("Mu")) %>% t() %>% as.data.frame()
Placenta <- OriExp70gene %>% select(.,contains("Pl")) %>% t() %>% as.data.frame()

TisCombine <- rbind(Brain,Liver,Muscle,Placenta)

TisCombine$Tissue <- rep(c("Brain","Liver","Muscle","Placenta"),each=12)

datPlot <- TisCombine[,c(3,6)]
colnames(datPlot) <- c("Gene","Tissue")

library(ggplot2)
library(ggsci)

ggplot(datPlot, aes(x=Tissue,y=Gene,fill=Tissue)) + geom_bar(stat = "summary", fun=mean, width = 0.5, color="black") + stat_summary(fun.data = 'mean_sdl', geom="errorbar", color="black", width=0.2, position = position_dodge(.9))+theme_bw() + scale_fill_aaas() + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("DLK1 expression (read count)")


##Changes of DLK1 expression in different periods
library(dplyr)

OriExp <- read.csv("F:\\Manuscript\\RNA\\Ori\\data\\Tissue_Stage_specific\\GeneExpressionAllSamples.fittedValue.csv", header = T, row.names = 1)

OriExp <- select(OriExp,contains("Br"))

gene <- c("ENSSSCG00000035805")

OriExpgene <- filter(OriExp, rownames(OriExp) %in% gene)

F40 <- OriExpgene %>% select(.,contains("40")) %>% t() %>% as.data.frame()
F70 <- OriExpgene %>% select(.,contains("70")) %>% t() %>% as.data.frame()
D1 <- OriExpgene %>% select(.,contains("115")) %>% t() %>% as.data.frame()
D168 <- OriExpgene %>% select(.,contains("168")) %>% t() %>% as.data.frame()

StgCombine <- rbind(F40,F70,D1,D168)

StgCombine$Stage <- factor(rep(c("F40","F70","D1","D168"),each=12),levels=c("F40","F70","D1","D168"))



datPlot <- StgCombine
colnames(datPlot) <- c("Gene","Stage")

library(ggplot2)
library(ggsci)

ggplot(datPlot, aes(x=Stage,y=Gene,fill=Stage)) + geom_bar(stat = "summary", fun=mean, width = 0.5, color="black") + stat_summary(fun.data = 'mean_sdl', geom="errorbar", color="black", width=0.2, position = position_dodge(.9))+theme_bw() + scale_fill_aaas() + theme(axis.title.x = element_blank(), legend.position = "none") + ylab("DLK1 expression (TPM)")

