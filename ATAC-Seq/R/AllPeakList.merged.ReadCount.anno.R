##针对Original genome peaks 的annotation

setwd("F:\\Manuscript\\ATAC\\data")

library(GenomicFeatures)
require(ChIPseeker)
library(clusterProfiler)

#txdb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl",host="www.ensembl.org")

txdb2 <- makeTxDbFromGFF("F:\\Resource\\Sus_scrofa.Sscrofa11.1.98.gff3")

peak_RC <- read.table("AllPeakList.merged.ReadCount.txt", header = T, row.names = 1)
peak <- readPeakFile("AllPeakList.merged.ReadCount.txt", as="GRanges")

peakAnno <- annotatePeak(peak, tssRegion = c(-3000,3000), TxDb = txdb2)

plotAnnoPie(peakAnno)

GeneID <- peakAnno@anno@elementMetadata@listData$geneId
GeneChr <- peakAnno@anno@elementMetadata@listData$geneChr
GeneStart <- peakAnno@anno@elementMetadata@listData$geneStart
GeneEnd <- peakAnno@anno@elementMetadata@listData$geneEnd
Distance <- peakAnno@anno@elementMetadata@listData$distance

##The feature types in the genome corresponding to peaks are obtained
Annotation <- peakAnno@detailGenomicAnnotation 

Peak_Gene <- data.frame(Peaks=rownames(peak_RC),gene=GeneID,Chr=GeneChr,Start=GeneStart,End=GeneEnd, DistanceToTSS=Distance)

Peak_Gene_Feature <- data.frame(Peak_Gene,Annotation)

write.csv(Peak_Gene_Feature,"AllPeakList.merged.ReadCount.PeakToGene_Feature.csv", quote = F, row.names = F)
