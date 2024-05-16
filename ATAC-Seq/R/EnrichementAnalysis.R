setwd("F:\\ATAC\\processFiles\\RNA-tissue_specific0711")
library(dplyr)
library(edgeR)
library(ggplot2)
library(stringr)

##gene expression matrix
dataExp <- read.csv("F:\\RNA\\DEG_Tissues_Stages\\GeneExpressionAllSamples.csv", header = T, row.names = 1)

##A list of genes involved in gene expression measurement in each tissue was extracted
group <- select(dataExp,contains("Pl"))
y <- DGEList(counts=group,genes=rownames(group))
#brain, liver, muscle
isexpr <- rowSums(cpm(y)>1) >= 15
#placenta
isexpr <- rowSums(cpm(y)>1) >= 10
y <- y[isexpr, , keep.lib.sizes=FALSE]
geneList <- y$genes
write.csv(geneList, "GeneList.Pl.csv", quote = F, row.names = TRUE)


##Find the gene that has the closest open peaks specific to each tissue
##in Linux
cd /mnt/scratch/quanjian/ATAC/peakCombine
join <(sort -k1,1 TisSpe.Pl.txt) <(sort -k1,1 peakList.bed)|awk 'OFS="\t"{print $2,$3,$4,$1}'|sort -k1,1n -k2,2n > TisSpe.Pl.bed
bedtools closest -a <(sort -k1,1 -k2,2n TisSpe.Pl.bed) -b <(sort -k1,1 -k2,2n ~/resources/RefGeneList.bed)|cut -f8|sort|uniq > TisSpe.Pl.closestGene.txt

join <(sort -k1,1 TisSpe.Br.txt) <(sort -k1,1 peakList.bed)|awk 'OFS="\t"{print $2,$3,$4,$1}'|sort -k1,1n -k2,2n > TisSpe.Br.bed
 bedtools closest -a <(sort -k1,1 -k2,2n TisSpe.Br.bed) -b <(sort -k1,1 -k2,2n ~/resources/RefGeneList.bed)|cut -f8|sort|uniq > TisSpe.Br.closestGene.txt

join <(sort -k1,1 TisSpe.Mu.txt) <(sort -k1,1 peakList.bed)|awk 'OFS="\t"{print $2,$3,$4,$1}'|sort -k1,1n -k2,2n > TisSpe.Mu.bed
bedtools closest -a <(sort -k1,1 -k2,2n TisSpe.Mu.bed) -b <(sort -k1,1 -k2,2n ~/resources/RefGeneList.bed)|cut -f8|sort|uniq > TisSpe.Mu.closestGene.txt

join <(sort -k1,1 TisSpe.Li.txt) <(sort -k1,1 peakList.bed)|awk 'OFS="\t"{print $2,$3,$4,$1}'|sort -k1,1n -k2,2n > TisSpe.Li.bed
bedtools closest -a <(sort -k1,1 -k2,2n TisSpe.Li.bed) -b <(sort -k1,1 -k2,2n ~/resources/RefGeneList.bed)|cut -f8|sort|uniq > TisSpe.Li.closestGene.txt


##enrichment analysis
setwd("F:\\Manuscript\\ATAC\\figure\\ori_enrichment")

geneListBr <- read.csv("GeneList.Br.csv", header = T, row.names = 1)
geneListLi <- read.csv("GeneList.Li.csv", header = T, row.names = 1)
geneListMu <- read.csv("GeneList.Mu.csv", header = T, row.names = 1)
geneListPl <- read.csv("GeneList.Pl.csv", header = T, row.names = 1)

BrainTis <- read.table("TisSpe.Br.closestGene.txt", header = F)
LiverTis <- read.table("TisSpe.Li.closestGene.txt", header = F)
MuscleTis <- read.table("TisSpe.Mu.closestGene.txt", header = F)
PlacentaTis <- read.table("TisSpe.Pl.closestGene.txt", header = F)

##Brain
brainSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListBr$genes, size = 388, replace = FALSE)
  overlap <- intersect(BrainTis$V1, tmp)
  brainSample <- c(brainSample, length(overlap))
}

##Liver
liverSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListLi$genes, size = 180, replace = FALSE)
  overlap <- intersect(LiverTis$V1, tmp)
  liverSample <- c(liverSample, length(overlap))
}

##Muscle
muscleSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListMu$genes, size = 131, replace = FALSE)
  overlap <- intersect(MuscleTis$V1, tmp)
  muscleSample <- c(muscleSample, length(overlap))
}

#Placenta
placentaSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListPl$genes, size = 268, replace = FALSE)
  overlap <- intersect(PlacentaTis$V1, tmp)
  placentaSample <- c(placentaSample, length(overlap))
}

##mean overlap of each tissue
enri_br <- mean(brainSample)
enri_li <- mean(liverSample)
enri_mu <- mean(muscleSample)
enri_pl <- mean(placentaSample)

##enrichment score
155/enri_br
11/enri_li
36/enri_mu
174/enri_pl

##pvalue 
pvalue_br <- sum(brainSample>155)/10000+1
pvalue_li <- sum(liverSample>11)/10000+1
pvalue_mu <- sum(muscleSample>36)/10000+1
pvalue_pl <- sum(placentaSample>174)/10000+1



####based on the peak annotation results
setwd("F:\\Manuscript\\ATAC\\figure\\ori_enrichment")

geneListBr <- read.csv("GeneList.Br.csv", header = T, row.names = 1)
geneListLi <- read.csv("GeneList.Li.csv", header = T, row.names = 1)
geneListMu <- read.csv("GeneList.Mu.csv", header = T, row.names = 1)
geneListPl <- read.csv("GeneList.Pl.csv", header = T, row.names = 1)

BrainTis <- read.table("TisSpe.Br.AnnotatedGene.txt", header = T)
LiverTis <- read.table("TisSpe.Li.AnnotatedGene.txt", header = T)
MuscleTis <- read.table("TisSpe.Mu.AnnotatedGene.txt", header = T)
PlacentaTis <- read.table("TisSpe.Pl.AnnotatedGene.txt", header = T)

##Brain
brainSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListBr$genes, size = 388, replace = FALSE)
  overlap <- intersect(BrainTis$., tmp)
  brainSample <- c(brainSample, length(overlap))
}

##Liver
liverSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListLi$genes, size = 180, replace = FALSE)
  overlap <- intersect(LiverTis$., tmp)
  liverSample <- c(liverSample, length(overlap))
}

##Muscle
muscleSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListMu$genes, size = 131, replace = FALSE)
  overlap <- intersect(MuscleTis$., tmp)
  muscleSample <- c(muscleSample, length(overlap))
}

#Placenta
placentaSample <- c()

for(i in 1:10000){
  tmp <- sample(geneListPl$genes, size = 268, replace = FALSE)
  overlap <- intersect(PlacentaTis$., tmp)
  placentaSample <- c(placentaSample, length(overlap))
}

##mean overlap of each tissue
enri_br <- mean(brainSample)
enri_li <- mean(liverSample)
enri_mu <- mean(muscleSample)
enri_pl <- mean(placentaSample)

##enrichment score
149/enri_br
15/enri_li
35/enri_mu
100/enri_pl

##pvalue 
pvalue_br <- (sum(brainSample>149)+1)/(10000+1)
pvalue_li <- (sum(liverSample>15)+1)/(10000+1)
pvalue_mu <- (sum(muscleSample>35)+1)/(10000+1)
pvalue_pl <- (sum(placentaSample>100)+1)/(10000+1)
