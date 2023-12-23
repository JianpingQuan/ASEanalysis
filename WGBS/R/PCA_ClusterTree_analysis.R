#cd /mnt/scratch/quanjian/BS/OriginalRef/MethylationExt
#module purge
#module load GCC/11.2.0  OpenMPI/4.1.1 R/4.1.2
#R

library(methylKit)
library(stringr)
library(dplyr)

files <- list.files(pattern="_OriRef.deduplicated.sortbyName.CpG_report.txt.gz")

file.list <- list()
for(i in 1:47){
  file.list[i] <- files[i]
}

#samples <- str_sub(files,1,-10)
samples <- str_sub(files,1,-50)
sample.list <- list()

for(i in 1:47){
  sample.list[i] <- samples[i]
}

samplename = samples

##compare between two group
treatment <- as.integer(grepl("_B", samplename))

###Multigroup comparison
treatment <- c()
for(i in 1:47){
 if(grepl("_B", samplename[i])){
  treatment[i] <- 0
 }else if(grepl("_L", samplename[i])){
  treatment[i] <- 1
 }else if(grepl("_M", samplename[i])){
  treatment[i] <- 2
 }else{
  treatment[i] <- 3
 }
}


myobj=methRead(file.list, sample.id=sample.list, assembly="sScr11", treatment=treatment, context="CpG",header=FALSE,pipeline="bismarkCytosineReport")

methylBase.obj <- unite(myobj, mc.cores = 10)

#View sample correlation
cor <- getCorrelation(methylBase.obj, method="pearson", plot=FALSE)

#samples clustering
pdf("MethylKit_Plot_correlation.pdf",height=4, width=8)
clusterSamples(methylBase.obj, dist="correlation", method="ward", plot=TRUE)
dev.off()

#PCA
pdf("MethylKit_Plot_PCA.pdf", heigth=6, width=6)
PCASamples(methylBase.obj)
dev.off()

##get result of prcomp
pca <- PCASamples(methylBase.obj, obj.return=TRUE)

##extract each PCs
m_df <- as.data.frame(pca$x)
summ <- summary(pca)
##get the variance importance of each PCs
pc <- as.data.frame(summ$importance)

xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

write.table(m_df, "SamplesPCA_prcomp_x.txt", quote = F, sep="\t")

write.table(pc, "SamplesPCA_variance_importance.txt", quote = F, sep="\t")