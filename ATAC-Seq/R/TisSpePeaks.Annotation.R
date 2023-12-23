library(dplyr)

dat <- read.csv("F:\\Manuscript\\ATAC\\figure\\ori_Tissue_stage_specific\\ReadCount\\ATAC-Tissue_specific.peaks.csv", header = T) %>% arrange(.,Var1)

AnnotationFile <- read.csv("F:\\Manuscript\\ATAC\\data\\AllPeakList.merged.ReadCount.PeakToGene_Feature.csv", header = T)

AnnotationFile_f <- filter(AnnotationFile, Peaks %in% dat$Var1) %>% arrange(.,Peaks)
dat$Gene <- AnnotationFile_f$gene


Br <- filter(dat, Tissue=="Br")[,4]
Li <- filter(dat, Tissue=="Li")[,4]
Mu <- filter(dat, Tissue=="Mu")[,4]
Pl <- filter(dat, Tissue=="Pl")[,4]


TisExpGene <- read.csv("F:\\ATAC\\processFiles\\RNA-tissue_specific0711\\RNA-Tissue_specific.genes.csv", header = T)


TisExpGene_Br <- filter(TisExpGene, Tissue=="Br")[,1]
TisExpGene_Li <- filter(TisExpGene, Tissue=="Li")[,1]
TisExpGene_Mu <- filter(TisExpGene, Tissue=="Mu")[,1]
TisExpGene_Pl <- filter(TisExpGene, Tissue=="Pl")[,1]


intersect(Br, TisExpGene_Br) %>% length()
intersect(Li, TisExpGene_Li) %>% length()
intersect(Mu, TisExpGene_Mu) %>% length()
intersect(Pl, TisExpGene_Pl) %>% length()



Br_uniq <- sort(Br) %>% unique() %>% as.data.frame()
Li_uniq <- sort(Li) %>% unique() %>% as.data.frame()
Mu_uniq <- sort(Mu) %>% unique() %>% as.data.frame()
Pl_uniq <- sort(Pl) %>% unique() %>% as.data.frame()


write.table(Br_uniq,"TisSpe.Br.AnnotatedGene.txt", quote = F, row.names = F)
write.table(Li_uniq,"TisSpe.Li.AnnotatedGene.txt", quote = F, row.names = F)
write.table(Mu_uniq,"TisSpe.Mu.AnnotatedGene.txt", quote = F, row.names = F)
write.table(Pl_uniq,"TisSpe.Pl.AnnotatedGene.txt", quote = F, row.names = F)

